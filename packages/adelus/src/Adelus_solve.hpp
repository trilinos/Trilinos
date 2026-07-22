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

#ifndef __ADELUS_SOLVE_HPP__
#define __ADELUS_SOLVE_HPP__

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "Adelus_defines.h"
#include "Adelus_macros.h"
#include "Adelus_mytime.hpp"
#include "Kokkos_Core.hpp"
#include "KokkosBlas3_gemm.hpp"

#define IBM_MPI_WRKAROUND2

#define SOSTATUSINT 32768

// Message tags
#define SOCOLTYPE (1<<14)
#define SOROWTYPE (1<<15)
#define SOHSTYPE (SOCOLTYPE + SOROWTYPE)

namespace Adelus {

//  Customized elimination on the rhs that I own
template<class ZView, class RHSView, class DView>
void elimination_rhs(int N, ZView& ptr2, RHSView& ptr3, DView& ptr4, int act_col) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
  Kokkos::parallel_for(Kokkos::RangePolicy<typename ZView::device_type::execution_space>(0,N), KOKKOS_LAMBDA (const int i) {
    ptr4(0,i) = ptr3(i)/ptr2(act_col);
    ptr3(i)   = ptr4(0,i);
  });
#else//OpenMP //Note: FOR loop is faster than parallel_for when execution_space is OpenMP
  for(int i=0;i<N;i++){
    ptr4(0,i) = ptr3(i)/ptr2(act_col);
    ptr3(i)   = ptr4(0,i);
  }
#endif
}

template<class HandleType, class ZViewType, class RHSViewType>
inline
void back_solve_rhs_pipelined_comm(HandleType& ahandle, ZViewType& Z, RHSViewType& RHS)
{
  using value_type      = typename ZViewType::value_type;
#ifdef PRINT_STATUS
  using execution_space = typename ZViewType::device_type::execution_space;
#endif
  using memory_space    = typename ZViewType::device_type::memory_space;
  using View2DType      = Kokkos::View<value_type**, Kokkos::LayoutLeft, memory_space>;

#if defined(ADELUS_HOST_PINNED_MEM_MPI) || defined(IBM_MPI_WRKAROUND2)
#if defined(KOKKOS_ENABLE_CUDA)
  using View2DHostPinnType = Kokkos::View<value_type**, Kokkos::LayoutLeft, Kokkos::CudaHostPinnedSpace>;//CudaHostPinnedSpace
#elif defined(KOKKOS_ENABLE_HIP)
  using View2DHostPinnType = Kokkos::View<value_type**, Kokkos::LayoutLeft, Kokkos::HIPHostPinnedSpace>;//HIPHostPinnedSpace
#endif
#endif

  MPI_Comm comm     = ahandle.get_comm();
  MPI_Comm col_comm = ahandle.get_col_comm();
  int me            = ahandle.get_myrank();
  int nprocs_row    = ahandle.get_nprocs_row();
  int nprocs_col    = ahandle.get_nprocs_col();
  int ncols_matrix  = ahandle.get_ncols_matrix();
  int my_rows       = ahandle.get_my_rows();
  int my_cols       = ahandle.get_my_cols();
  int my_first_row  = ahandle.get_my_first_row();
  int my_first_col  = ahandle.get_my_first_col();
  int nrhs          = ahandle.get_nrhs();
  int my_rhs        = ahandle.get_my_rhs();

  int  j;         // loop counters
  int end_row;    // row num to end column operations
  int bytes[16];  // number of bytes in messages
  int root;       // root processor for fanout
  int type[16];   // mesage type for messages
  int dest[16];   // dest for message sends

  int one = 1;

  value_type d_one = 1.0;
  value_type d_min_one = -1.0;

  int j2;

  int n_rhs_this; // num rhs that I currently own
  int col_offset; // which processor starts the pipeline
  int my_pos;     // my position in the new linup
  int extra;      // extra loop to realign data after pipeline
  int act_col;    // act this column (that I own)
  int on_col;     // on this collection of rhs's
  int global_col; // global col number for act_col
  int max_bytes;  // max number of bytes of rhs I can receive

  int my_col_id, my_row_id, id_temp;
  int dest_right, dest_left;

#ifdef GET_TIMING
  double t1,t2;
  double allocviewtime,eliminaterhstime,bcastrowtime,updrhstime,xchgrhstime;
  double totalsolvetime;
#if defined(ADELUS_HOST_PINNED_MEM_MPI) && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
  double copyhostpinnedtime;
#endif
#endif

  MPI_Request msgrequest;
  MPI_Status msgstatus;

#ifdef GET_TIMING
  t2 = MPI_Wtime();
#endif

  // find left, right destination procs

  my_col_id = mesh_col(me);
  my_row_id = mesh_row(me);

  id_temp = my_col_id + 1;
  if (id_temp >= nprocs_row) id_temp = 0;
  dest_right = proc_num(my_row_id,id_temp);

  id_temp = my_col_id - 1;
  if (id_temp < 0) id_temp = nprocs_row-1;
  dest_left = proc_num(my_row_id,id_temp);

  // set j2 to be first column in last group of columns
  max_bytes = nrhs/nprocs_row;
  if (nrhs%nprocs_row > 0) max_bytes++;
  max_bytes = max_bytes*sizeof(ADELUS_DATA_TYPE)*my_rows;

#ifdef GET_TIMING
  allocviewtime=eliminaterhstime=bcastrowtime=updrhstime=xchgrhstime=0.0;
#if defined(ADELUS_HOST_PINNED_MEM_MPI) && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
  copyhostpinnedtime=0.0;
#endif

  t1 = MPI_Wtime();
#endif

  View2DType row1( "row1", one, nrhs );   // row1: diagonal row (temp variables)
#if (defined(ADELUS_HOST_PINNED_MEM_MPI) || defined(IBM_MPI_WRKAROUND2)) && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
  View2DHostPinnType h_row2( "h_row2", my_rows, max_bytes/sizeof(ADELUS_DATA_TYPE)/my_rows );
#else
  View2DType row2( "row2", my_rows, max_bytes/sizeof(ADELUS_DATA_TYPE)/my_rows );
#endif
#if (defined(ADELUS_HOST_PINNED_MEM_MPI) || defined(IBM_MPI_WRKAROUND2)) && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
  View2DHostPinnType h_row1( "h_row1", one, nrhs );
  View2DHostPinnType h_rhs ( "h_rhs",  my_rows, nrhs );
#endif

  Kokkos::fence();

#ifdef GET_TIMING
  allocviewtime += (MPI_Wtime()-t1);
#endif

  n_rhs_this = my_rhs;// why "n_rhs_this" is changing after the first iteration (need to as Joe???)
  j2 = ncols_matrix-1;
  col_offset = (j2%nprocs_row);
  my_pos = my_first_col - col_offset;
  if (my_pos < 0) my_pos += nprocs_row;
  extra = (nprocs_row - (col_offset-1))%nprocs_row;

  act_col = my_cols-1;
  if (my_pos != 0) act_col++;

  on_col = my_pos;

#ifdef PRINT_STATUS
  printf("Rank %i -- back_solve6() Begin back solve, execution_space %s, memory_space %s\n",me, typeid(execution_space).name(), typeid(memory_space).name());
#endif

  for (j = j2; j >= 1-nprocs_row-extra; j--) {

    if ((j+nprocs_row-1 >= 0) && (n_rhs_this > 0)) {

      if ((act_col < my_cols) && (act_col >= 0)) {

        global_col = act_col*nprocs_row + my_first_col;

        end_row = global_col/nprocs_col;
        if (my_first_row <= global_col%nprocs_col) ++end_row;

        // do an elimination step on the rhs that I own

        //auto ptr2_view = subview(Z, end_row-1, Kokkos::ALL());

        root = row_owner(global_col);

        if (me == root) {
#ifdef GET_TIMING
          t1 = MPI_Wtime();
#endif
          auto ptr2_view = subview(Z,   end_row-1, Kokkos::ALL());
          auto ptr3_view = subview(RHS, end_row-1, Kokkos::make_pair(0, n_rhs_this));
          elimination_rhs(n_rhs_this, ptr2_view, ptr3_view, row1, act_col);//note: row1 = ptr4
//Workaround for GPU-aware MPI issue on HIP: fence here to make sure "elimination_rhs" completes before MPI_Bcast
#if defined(KOKKOS_ENABLE_HIP) && !defined(ADELUS_HOST_PINNED_MEM_MPI)
          Kokkos::fence();
#endif
          end_row--;
#ifdef GET_TIMING
          eliminaterhstime += (MPI_Wtime()-t1);
#endif
        }

#if defined(ADELUS_HOST_PINNED_MEM_MPI) && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
#ifdef GET_TIMING
        t1 = MPI_Wtime();
#endif
        Kokkos::deep_copy(h_row1,row1);
#ifdef GET_TIMING
        copyhostpinnedtime += (MPI_Wtime()-t1);
#endif
#endif

#ifdef GET_TIMING
        t1 = MPI_Wtime();
#endif
        bytes[0] = n_rhs_this*sizeof(ADELUS_DATA_TYPE);
        type[0]  = SOCOLTYPE+j;

        //MPI_Bcast((char *) row1, bytes[0], MPI_CHAR, mesh_row(root), col_comm);
#if defined(ADELUS_HOST_PINNED_MEM_MPI) && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
        MPI_Bcast(reinterpret_cast<char *>(h_row1.data()), bytes[0], MPI_CHAR, mesh_row(root), col_comm);
#else //GPU-aware MPI
        MPI_Bcast(reinterpret_cast<char *>(row1.data()), bytes[0], MPI_CHAR, mesh_row(root), col_comm);
#endif
        // added this barrier for CPLANT operation

        MPI_Barrier(col_comm);
#ifdef GET_TIMING
        bcastrowtime += (MPI_Wtime()-t1);
#endif

#if defined(ADELUS_HOST_PINNED_MEM_MPI) && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
#ifdef GET_TIMING
        t1 = MPI_Wtime();
#endif
        Kokkos::deep_copy(row1,h_row1);
#ifdef GET_TIMING
        copyhostpinnedtime += (MPI_Wtime()-t1);
#endif
#endif

#ifdef GET_TIMING
        t1 = MPI_Wtime();
#endif

        auto A_view = subview(Z,    Kokkos::make_pair(0, end_row), Kokkos::make_pair(act_col, act_col+one));
        auto C_view = subview(RHS,  Kokkos::make_pair(0, end_row), Kokkos::make_pair(0, n_rhs_this));
        auto B_view = subview(row1, Kokkos::ALL(), Kokkos::make_pair(0, n_rhs_this));

        KokkosBlas::gemm("N","N",d_min_one,
                         A_view,
                         B_view,
                         d_one,
                         C_view);
//Workaround for GPU-aware MPI issue on HIP: fence here to make sure "gemm" completes before MPI_Send
#if defined(KOKKOS_ENABLE_HIP) && !defined(ADELUS_HOST_PINNED_MEM_MPI)
        Kokkos::fence();
#endif
#ifdef GET_TIMING
        updrhstime += (MPI_Wtime()-t1);
#endif
      }
    }

#ifdef GET_TIMING
    t1 = MPI_Wtime();
#endif
    if (j != 1-nprocs_row-extra) {
      dest[0] = dest_right;
      if (me != dest[0]) {
        bytes[0] = max_bytes;
        type[0]  = SOROWTYPE+j;

#if (defined(ADELUS_HOST_PINNED_MEM_MPI) || defined(IBM_MPI_WRKAROUND2)) && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
        MPI_Irecv(reinterpret_cast<char *>(h_row2.data()), bytes[0], MPI_CHAR, MPI_ANY_SOURCE, type[0], comm, &msgrequest);
#else
        MPI_Irecv(reinterpret_cast<char *>(  row2.data()), bytes[0], MPI_CHAR, MPI_ANY_SOURCE, type[0], comm, &msgrequest);
#endif

        n_rhs_this = bytes[0]/sizeof(ADELUS_DATA_TYPE)/my_rows;

#if (defined(ADELUS_HOST_PINNED_MEM_MPI) || defined(IBM_MPI_WRKAROUND2)) && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
        Kokkos::deep_copy(subview(h_rhs, Kokkos::ALL(), Kokkos::make_pair(0, n_rhs_this)),
                          subview(RHS,   Kokkos::ALL(), Kokkos::make_pair(0, n_rhs_this)));
#endif

        dest[1]  = dest_left;
        bytes[1] = n_rhs_this * sizeof(ADELUS_DATA_TYPE) * my_rows;
        type[1]  = SOROWTYPE+j;

#if (defined(ADELUS_HOST_PINNED_MEM_MPI) || defined(IBM_MPI_WRKAROUND2)) && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
        MPI_Send(reinterpret_cast<char *>(h_rhs.data()), bytes[1], MPI_CHAR, dest[1], type[1], comm);
#else //GPU-aware MPI
        MPI_Send(reinterpret_cast<char *>(RHS.data()), bytes[1], MPI_CHAR, dest[1], type[1], comm);
#endif

        MPI_Wait(&msgrequest,&msgstatus);

        // Copy row2 -> rhs
        int blas_length = n_rhs_this*my_rows;
#if (defined(ADELUS_HOST_PINNED_MEM_MPI) || defined(IBM_MPI_WRKAROUND2)) && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)) //Use memcpy for now, can use deep_copy in the future //deep_copy is slower than BLAS XCOPY
#if defined(KOKKOS_ENABLE_CUDA)
        //Kokkos::deep_copy(subview(RHS, Kokkos::ALL(), Kokkos::make_pair(0, n_rhs_this)), subview(h_row2, Kokkos::ALL(), Kokkos::make_pair(0, n_rhs_this)));
        cudaMemcpy(reinterpret_cast<ADELUS_DATA_TYPE *>(RHS.data()), reinterpret_cast<ADELUS_DATA_TYPE *>(h_row2.data()), blas_length*sizeof(ADELUS_DATA_TYPE), cudaMemcpyHostToDevice);
#elif defined(KOKKOS_ENABLE_HIP)
        hipMemcpy(reinterpret_cast<ADELUS_DATA_TYPE *>(RHS.data()), reinterpret_cast<ADELUS_DATA_TYPE *>(h_row2.data()), blas_length*sizeof(ADELUS_DATA_TYPE), hipMemcpyHostToDevice);
#endif
#else
#if defined(KOKKOS_ENABLE_CUDA)
        cudaMemcpy(reinterpret_cast<ADELUS_DATA_TYPE *>(RHS.data()), reinterpret_cast<ADELUS_DATA_TYPE *>(row2.data()), blas_length*sizeof(ADELUS_DATA_TYPE), cudaMemcpyDeviceToDevice);
#elif defined(KOKKOS_ENABLE_HIP)
        hipMemcpy(reinterpret_cast<ADELUS_DATA_TYPE *>(RHS.data()), reinterpret_cast<ADELUS_DATA_TYPE *>(row2.data()), blas_length*sizeof(ADELUS_DATA_TYPE), hipMemcpyDeviceToDevice);
#else
        memcpy(reinterpret_cast<ADELUS_DATA_TYPE *>(RHS.data()), reinterpret_cast<ADELUS_DATA_TYPE *>(row2.data()), blas_length*sizeof(ADELUS_DATA_TYPE));
#endif
#endif
      }
      on_col++;
      if (on_col >= nprocs_row) {
        on_col = 0;
        act_col--;
      }
    }
#ifdef GET_TIMING
    xchgrhstime += (MPI_Wtime()-t1);
#endif

  }

#ifdef GET_TIMING
  totalsolvetime = MPI_Wtime() - t2;
#endif
#ifdef GET_TIMING
  showtime(ahandle.get_comm_id(), comm, me, ahandle.get_nprocs_cube(), "Time to alloc view", &allocviewtime);
  showtime(ahandle.get_comm_id(), comm, me, ahandle.get_nprocs_cube(), "Time to eliminate rhs",&eliminaterhstime);
  showtime(ahandle.get_comm_id(), comm, me, ahandle.get_nprocs_cube(), "Time to bcast temp row",&bcastrowtime);
  showtime(ahandle.get_comm_id(), comm, me, ahandle.get_nprocs_cube(), "Time to update rhs",&updrhstime);
#if defined(ADELUS_HOST_PINNED_MEM_MPI) && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
  showtime(ahandle.get_comm_id(), comm, me, ahandle.get_nprocs_cube(), "Time to copy host pinned mem <--> dev mem",&copyhostpinnedtime);
#endif
  showtime(ahandle.get_comm_id(), comm, me, ahandle.get_nprocs_cube(), "Time to xchg rhs",&xchgrhstime);
  showtime(ahandle.get_comm_id(), comm, me, ahandle.get_nprocs_cube(), "Total time in solve",&totalsolvetime);
#endif
}

template<class HandleType, class ZViewType, class RHSViewType>
inline
void back_solve_currcol_bcast(HandleType& ahandle, ZViewType& Z, RHSViewType& RHS)
{
  using value_type      = typename ZViewType::value_type;
#ifdef PRINT_STATUS
  using execution_space = typename ZViewType::device_type::execution_space;
#endif
  using memory_space    = typename ZViewType::device_type::memory_space;
  using View2DType      = Kokkos::View<value_type**, Kokkos::LayoutLeft, memory_space>;

#if defined(ADELUS_HOST_PINNED_MEM_MPI)
#if defined(KOKKOS_ENABLE_CUDA)
  using View2DHostPinnType = Kokkos::View<value_type**, Kokkos::LayoutLeft, Kokkos::CudaHostPinnedSpace>;//CudaHostPinnedSpace
#elif defined(KOKKOS_ENABLE_HIP)
  using View2DHostPinnType = Kokkos::View<value_type**, Kokkos::LayoutLeft, Kokkos::HIPHostPinnedSpace>;//HIPHostPinnedSpace
#endif
#endif

#if defined(GET_TIMING) || defined(PRINT_STATUS)
  int me            = ahandle.get_myrank();
#ifdef GET_TIMING
  MPI_Comm comm     = ahandle.get_comm();
#endif
#endif
  MPI_Comm col_comm = ahandle.get_col_comm();
  MPI_Comm row_comm = ahandle.get_row_comm();
  int myrow         = ahandle.get_myrow();
  int mycol         = ahandle.get_mycol();
  int nprocs_row    = ahandle.get_nprocs_row();
  int nprocs_col    = ahandle.get_nprocs_col();
  int ncols_matrix  = ahandle.get_ncols_matrix();
  int my_rows       = ahandle.get_my_rows();
  int my_rhs        = ahandle.get_my_rhs();

  value_type d_one = 1.0;
  value_type d_min_one = -1.0;

#ifdef GET_TIMING
  double t1,t2;
  double allocviewtime,eliminaterhstime,bcastrowtime,updrhstime,bcastcoltime,copycoltime;
  double totalsolvetime;
#if defined(ADELUS_HOST_PINNED_MEM_MPI) && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
  double copyhostpinnedtime;
#endif
#endif

#ifdef GET_TIMING
  t2 = MPI_Wtime();
#endif

#ifdef GET_TIMING
  allocviewtime=eliminaterhstime=bcastrowtime=updrhstime=bcastcoltime=copycoltime=0.0;
#if defined(ADELUS_HOST_PINNED_MEM_MPI) && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
  copyhostpinnedtime=0.0;
#endif

  t1 = MPI_Wtime();
#endif

  View2DType curr_col( "curr_col", my_rows, 1 ); //current column
  View2DType rhs_row ( "rhs_row", 1, my_rhs );   //current row of RHS to hold the elimination results (i.e row of solution)
#if defined(ADELUS_HOST_PINNED_MEM_MPI) && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
  View2DHostPinnType h_curr_col( "h_curr_col", my_rows, 1 );
  View2DHostPinnType h_rhs_row( "h_rhs_row", 1, my_rhs );
#endif

  //Kokkos::fence();//NOTE: Should we need this?

#ifdef GET_TIMING
  allocviewtime += (MPI_Wtime()-t1);
#endif

#ifdef PRINT_STATUS
  printf("Rank %i -- back_solve6() Begin back solve, execution_space %s, memory_space %s\n",me, typeid(execution_space).name(), typeid(memory_space).name());
#endif

  for (int k = ncols_matrix-1; k >= 0; k--) {
    int k_row = k%nprocs_col;//proc. id (in the col_comm) having global k
    int k_col = k%nprocs_row;//proc. id (in the row_comm) having global k
    int end_row = k/nprocs_col;
    if (myrow <= k_row) end_row++;

#ifdef GET_TIMING
    t1 = MPI_Wtime();
#endif
    //Step 1: copy the current column of Z to a temporary view
    if (mycol == k_col) { //only deep_copy if holding the current column
      Kokkos::deep_copy( Kokkos::subview(curr_col, Kokkos::make_pair(0, end_row), 0),
                         Kokkos::subview(Z, Kokkos::make_pair(0, end_row), k/nprocs_row) );
    }
#ifdef GET_TIMING
    copycoltime += (MPI_Wtime()-t1);
#endif

#if defined(ADELUS_HOST_PINNED_MEM_MPI) && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
#ifdef GET_TIMING
    t1 = MPI_Wtime();
#endif
    Kokkos::deep_copy(h_curr_col,curr_col);
#ifdef GET_TIMING
    copyhostpinnedtime += (MPI_Wtime()-t1);
#endif
#endif

#ifdef GET_TIMING
    t1 = MPI_Wtime();
#endif
    //Step 2: broadcast the current column to all ranks in the row_comm
#if defined(ADELUS_HOST_PINNED_MEM_MPI) && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
    MPI_Bcast(reinterpret_cast<char *>(h_curr_col.data()), end_row*sizeof(ADELUS_DATA_TYPE), MPI_CHAR, k_col, row_comm);
#else //GPU-aware MPI
    MPI_Bcast(reinterpret_cast<char *>(curr_col.data()), end_row*sizeof(ADELUS_DATA_TYPE), MPI_CHAR, k_col, row_comm);
#endif
#ifdef GET_TIMING
    bcastcoltime += (MPI_Wtime()-t1);
#endif

#if defined(ADELUS_HOST_PINNED_MEM_MPI) && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
#ifdef GET_TIMING
    t1 = MPI_Wtime();
#endif
    Kokkos::deep_copy(curr_col,h_curr_col);
#ifdef GET_TIMING
    copyhostpinnedtime += (MPI_Wtime()-t1);
#endif
#endif

#ifdef GET_TIMING
    t1 = MPI_Wtime();
#endif
    //Step 3: do rhs elimination to get solution x
    if (myrow == k_row) {//only on ranks having row k
      if (my_rhs > 0) {  //only on ranks having some rhs
        auto sub_curr_col = Kokkos::subview(curr_col, end_row-1, Kokkos::ALL());
        auto sub_rhs      = Kokkos::subview(RHS,      end_row-1, Kokkos::make_pair(0, my_rhs));
        int act_col = 0;
        elimination_rhs(my_rhs, sub_curr_col, sub_rhs, rhs_row, act_col); Kokkos::fence();
        end_row--;//do not count the eliminated row in Step 5
      }
    }
#ifdef GET_TIMING
    eliminaterhstime += (MPI_Wtime()-t1);
#endif

    //MPI_Barrier(comm);//NOTE: Should we need this?

    if (my_rhs > 0) { //only on ranks having rhs
      if (k >= 1) {//still have row(s) to do rhs updates with elimination results
#if defined(ADELUS_HOST_PINNED_MEM_MPI) && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
#ifdef GET_TIMING
        t1 = MPI_Wtime();
#endif
        Kokkos::deep_copy(h_rhs_row,rhs_row);
#ifdef GET_TIMING
        copyhostpinnedtime += (MPI_Wtime()-t1);
#endif
#endif

#ifdef GET_TIMING
        t1 = MPI_Wtime();
#endif
        //Step 4: broadcast elimination results to all ranks in col_comm
#if defined(ADELUS_HOST_PINNED_MEM_MPI) && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
        MPI_Bcast(reinterpret_cast<char *>(h_rhs_row.data()), my_rhs*sizeof(ADELUS_DATA_TYPE), MPI_CHAR, k_row, col_comm);
#else //GPU-aware MPI
        MPI_Bcast(reinterpret_cast<char *>(rhs_row.data()), my_rhs*sizeof(ADELUS_DATA_TYPE), MPI_CHAR, k_row, col_comm);
#endif

#ifdef GET_TIMING
        bcastrowtime += (MPI_Wtime()-t1);
#endif

#if defined(ADELUS_HOST_PINNED_MEM_MPI) && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
#ifdef GET_TIMING
        t1 = MPI_Wtime();
#endif
        Kokkos::deep_copy(rhs_row,h_rhs_row);
#ifdef GET_TIMING
        copyhostpinnedtime += (MPI_Wtime()-t1);
#endif
#endif

#ifdef GET_TIMING
        t1 = MPI_Wtime();
#endif
        //Step 5: call gemm to update RHS with partial solution
        auto A_view = Kokkos::subview(curr_col, Kokkos::make_pair(0, end_row), Kokkos::ALL());
        auto B_view = Kokkos::subview(rhs_row,  Kokkos::ALL(), Kokkos::make_pair(0, my_rhs));
        auto C_view = Kokkos::subview(RHS,      Kokkos::make_pair(0, end_row), Kokkos::make_pair(0, my_rhs));

        KokkosBlas::gemm("N","N",d_min_one, A_view, B_view, d_one, C_view); Kokkos::fence();
#ifdef GET_TIMING
        updrhstime += (MPI_Wtime()-t1);
#endif
      }//end of (k >= 1)
    }//end of (my_rhs > 0)

    //MPI_Barrier(comm);//NOTE: Should we need this?

  }//end of for (int k = ncols_matrix-1; k >= 0; k--)

#ifdef GET_TIMING
  totalsolvetime = MPI_Wtime() - t2;
#endif
#ifdef GET_TIMING
  showtime(ahandle.get_comm_id(), comm, me, ahandle.get_nprocs_cube(), "Time to alloc view", &allocviewtime);
  showtime(ahandle.get_comm_id(), comm, me, ahandle.get_nprocs_cube(), "Time to copy matrix column",&copycoltime);
  showtime(ahandle.get_comm_id(), comm, me, ahandle.get_nprocs_cube(), "Time to bcast matrix column",&bcastcoltime);
  showtime(ahandle.get_comm_id(), comm, me, ahandle.get_nprocs_cube(), "Time to eliminate rhs",&eliminaterhstime);
  showtime(ahandle.get_comm_id(), comm, me, ahandle.get_nprocs_cube(), "Time to bcast temp row",&bcastrowtime);
  showtime(ahandle.get_comm_id(), comm, me, ahandle.get_nprocs_cube(), "Time to update rhs",&updrhstime);
#if defined(ADELUS_HOST_PINNED_MEM_MPI) && (defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP))
  showtime(ahandle.get_comm_id(), comm, me, ahandle.get_nprocs_cube(), "Time to copy host pinned mem <--> dev mem",&copyhostpinnedtime);
#endif
  showtime(ahandle.get_comm_id(), comm, me, ahandle.get_nprocs_cube(), "Total time in solve",&totalsolvetime);
#endif
}

template<class HandleType, class ZViewType, class RHSViewType>
inline
void back_solve6(HandleType& ahandle, ZViewType& Z, RHSViewType& RHS)
{
  //if (ahandle.get_nrhs() <= ahandle.get_nprocs_row()) {
  //  back_solve_rhs_pipelined_comm(ahandle, Z, RHS);
  //}
  //else {
    back_solve_currcol_bcast(ahandle, Z, RHS);
  //}
}

}//namespace Adelus

#endif
