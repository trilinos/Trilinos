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

#ifndef __ADELUS_SOLVE_HPP__
#define __ADELUS_SOLVE_HPP__

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "Adelus_defines.h"
#include "Adelus_macros.h"
#include "Adelus_pcomm.hpp"
#include "Adelus_mytime.hpp"
#include "Kokkos_Core.hpp"
#include "KokkosBlas3_gemm.hpp"

#define IBM_MPI_WRKAROUND2

extern int me;

extern int ncols_matrix;  // number of cols in the matrix

extern int nprocs_col;    // num of procs to which a col is assigned
extern int nprocs_row;    // num of procs to which a row is assigned

extern int my_first_col;  // proc position in a col
extern int my_first_row;  // proc position in a row

extern int my_rows;       // num of rows I own
extern int my_cols;       // num of cols I own

extern int nrhs;          // number of right hand sides
extern int my_rhs;        // number of right hand sides that I own

extern MPI_Comm col_comm;


#define SOSTATUSINT 32768

// Message tags
#define SOCOLTYPE (1<<14)
#define SOROWTYPE (1<<15)
#define SOHSTYPE (SOCOLTYPE + SOROWTYPE)

namespace Adelus {

//  Customized elimination on the rhs that I own	
template<class ZDView, class RView>
void elimination_rhs(int N, ZDView& ptr3, ZDView& ptr2, RView& ptr4, int act_col) {
#ifdef KOKKOS_ENABLE_CUDA
  Kokkos::parallel_for(Kokkos::RangePolicy<typename ZDView::device_type::execution_space>(0,N), KOKKOS_LAMBDA (const int i) {
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

template<class ZDView>
inline
void back_solve6(ZDView& ZV)
{
  typedef typename ZDView::value_type value_type;
#ifdef PRINT_STATUS
  typedef typename ZDView::device_type::execution_space execution_space;
#endif
  typedef typename ZDView::device_type::memory_space memory_space;
  typedef Kokkos::View<value_type**, Kokkos::LayoutLeft, memory_space> ViewMatrixType;

#if (defined(CUDA_HOST_PINNED_MPI) || defined(IBM_MPI_WRKAROUND2)) && defined(KOKKOS_ENABLE_CUDA)
  typedef Kokkos::View<value_type**, Kokkos::LayoutLeft, Kokkos::CudaHostPinnedSpace> View2DHostPinnType;//CudaHostPinnedSpace
#endif

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
//#if defined(CUDA_HOST_PINNED_MPI) && defined(KOKKOS_ENABLE_CUDA)
//  double copyhostpinnedtime;
//#endif
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
//#if defined(CUDA_HOST_PINNED_MPI) && defined(KOKKOS_ENABLE_CUDA)
//  copyhostpinnedtime=0.0;
//#endif

  t1 = MPI_Wtime();
#endif

  ViewMatrixType row1( "row1", one, nrhs );   // row1: diagonal row (temp variables)
#if (defined(CUDA_HOST_PINNED_MPI) || defined(IBM_MPI_WRKAROUND2)) && defined(KOKKOS_ENABLE_CUDA)
  View2DHostPinnType h_row2( "h_row2", my_rows, max_bytes/sizeof(ADELUS_DATA_TYPE)/my_rows );
#else
  ViewMatrixType row2( "row2", my_rows, max_bytes/sizeof(ADELUS_DATA_TYPE)/my_rows );
#endif
#if (defined(CUDA_HOST_PINNED_MPI) || defined(IBM_MPI_WRKAROUND2)) && defined(KOKKOS_ENABLE_CUDA)
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

        //auto ptr2_view = subview(ZV, end_row-1, Kokkos::ALL());

        root = row_owner(global_col);

        if (me == root) {
#ifdef GET_TIMING
          t1 = MPI_Wtime();
#endif
          auto ptr2_view = subview(ZV, end_row-1, Kokkos::ALL());
          auto ptr3_view = subview(ZV, end_row-1, Kokkos::make_pair(my_cols, my_cols+n_rhs_this));
          elimination_rhs(n_rhs_this, ptr3_view, ptr2_view, row1, act_col);//note: row1 = ptr4
          end_row--;
#ifdef GET_TIMING
          eliminaterhstime += (MPI_Wtime()-t1);
#endif
        }

//#if defined(CUDA_HOST_PINNED_MPI) && defined(KOKKOS_ENABLE_CUDA)
//#ifdef GET_TIMING
//        t1 = MPI_Wtime();
//#endif
//        Kokkos::deep_copy(h_row1,row1);
//#ifdef GET_TIMING
//        copyhostpinnedtime += (MPI_Wtime()-t1);
//#endif
//#endif

#ifdef GET_TIMING
        t1 = MPI_Wtime();
#endif
        bytes[0] = n_rhs_this*sizeof(ADELUS_DATA_TYPE);
        type[0]  = SOCOLTYPE+j;

        //MPI_Bcast((char *) row1, bytes[0], MPI_CHAR, mesh_row(root), col_comm);
//#if defined(CUDA_HOST_PINNED_MPI) && defined(KOKKOS_ENABLE_CUDA)
//        MPI_Bcast(reinterpret_cast<char *>(h_row1.data()), bytes[0], MPI_CHAR, mesh_row(root), col_comm);
//#else //CUDA-aware MPI -- Note: Looks like MPI_Bcast is still working well with device (cuda) pointers (and faster than using cuda host pinned memory)
        MPI_Bcast(reinterpret_cast<char *>(row1.data()), bytes[0], MPI_CHAR, mesh_row(root), col_comm);		
//#endif
        // added this barrier for CPLANT operation

        MPI_Barrier(col_comm);
#ifdef GET_TIMING
        bcastrowtime += (MPI_Wtime()-t1);
#endif

//#if defined(CUDA_HOST_PINNED_MPI) && defined(KOKKOS_ENABLE_CUDA)
//#ifdef GET_TIMING
//        t1 = MPI_Wtime();
//#endif
//        Kokkos::deep_copy(row1,h_row1);
//#ifdef GET_TIMING
//        copyhostpinnedtime += (MPI_Wtime()-t1);
//#endif
//#endif

#ifdef GET_TIMING
        t1 = MPI_Wtime();
#endif

        auto A_view = subview(ZV, Kokkos::make_pair(0, end_row), Kokkos::make_pair(act_col, act_col+one));
        auto C_view = subview(ZV, Kokkos::make_pair(0, end_row), Kokkos::make_pair(my_cols, my_cols+n_rhs_this));
        auto B_view = subview(row1, Kokkos::ALL(), Kokkos::make_pair(0, n_rhs_this));

        KokkosBlas::gemm("N","N",d_min_one,
                         A_view,
                         B_view,
                         d_one,
                         C_view);

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

#if (defined(CUDA_HOST_PINNED_MPI) || defined(IBM_MPI_WRKAROUND2)) && defined(KOKKOS_ENABLE_CUDA)
        MPI_Irecv(reinterpret_cast<char *>(h_row2.data()), bytes[0], MPI_CHAR, MPI_ANY_SOURCE, type[0], MPI_COMM_WORLD, &msgrequest);
#else
        MPI_Irecv(reinterpret_cast<char *>(  row2.data()), bytes[0], MPI_CHAR, MPI_ANY_SOURCE, type[0], MPI_COMM_WORLD, &msgrequest);
#endif

        n_rhs_this = bytes[0]/sizeof(ADELUS_DATA_TYPE)/my_rows;

#if (defined(CUDA_HOST_PINNED_MPI) || defined(IBM_MPI_WRKAROUND2)) && defined(KOKKOS_ENABLE_CUDA)
        Kokkos::deep_copy(subview(h_rhs, Kokkos::ALL(), Kokkos::make_pair(0, n_rhs_this)), subview(ZV, Kokkos::ALL(), Kokkos::make_pair(my_cols, my_cols+n_rhs_this)));
#endif

        dest[1]  = dest_left;
        bytes[1] = n_rhs_this * sizeof(ADELUS_DATA_TYPE) * my_rows;
        type[1]  = SOROWTYPE+j;

#if (defined(CUDA_HOST_PINNED_MPI) || defined(IBM_MPI_WRKAROUND2)) && defined(KOKKOS_ENABLE_CUDA)
        MPI_Send(reinterpret_cast<char *>(h_rhs.data()), bytes[1], MPI_CHAR, dest[1], type[1], MPI_COMM_WORLD);
#else //CUDA-aware MPI
        MPI_Send(reinterpret_cast<char *>(ZV.data()+my_rows*my_cols), bytes[1], MPI_CHAR, dest[1], type[1], MPI_COMM_WORLD);
#endif

        MPI_Wait(&msgrequest,&msgstatus);

        // Copy row2 -> rhs
        int blas_length = n_rhs_this*my_rows;
#if (defined(CUDA_HOST_PINNED_MPI) || defined(IBM_MPI_WRKAROUND2)) && defined(KOKKOS_ENABLE_CUDA)//Use memcpy for now, can use deep_copy in the future //deep_copy is slower than BLAS XCOPY
        //Kokkos::deep_copy(subview(ZV, Kokkos::ALL(), Kokkos::make_pair(my_cols, my_cols+n_rhs_this)), subview(h_row2, Kokkos::ALL(), Kokkos::make_pair(0, n_rhs_this)));
        cudaMemcpy(reinterpret_cast<ADELUS_DATA_TYPE *>(ZV.data()+my_rows*my_cols), reinterpret_cast<ADELUS_DATA_TYPE *>(h_row2.data()), blas_length*sizeof(ADELUS_DATA_TYPE), cudaMemcpyHostToDevice);
#else
#ifdef KOKKOS_ENABLE_CUDA
        cudaMemcpy(reinterpret_cast<ADELUS_DATA_TYPE *>(ZV.data()+my_rows*my_cols), reinterpret_cast<ADELUS_DATA_TYPE *>(row2.data()), blas_length*sizeof(ADELUS_DATA_TYPE), cudaMemcpyDeviceToDevice);
#else
        memcpy(reinterpret_cast<ADELUS_DATA_TYPE *>(ZV.data()+my_rows*my_cols), reinterpret_cast<ADELUS_DATA_TYPE *>(row2.data()), blas_length*sizeof(ADELUS_DATA_TYPE));
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
  showtime("Time to alloc view",&allocviewtime);
  showtime("Time to eliminate rhs",&eliminaterhstime);
  showtime("Time to bcast temp row",&bcastrowtime);
  showtime("Time to update rhs",&updrhstime);
//#if defined(CUDA_HOST_PINNED_MPI) && defined(KOKKOS_ENABLE_CUDA)
//  showtime("Time to copy host pinned mem <--> dev mem",&copyhostpinnedtime);   
//#endif
  showtime("Time to xchg rhs",&xchgrhstime);
  showtime("Total time in solve",&totalsolvetime);
#endif
}

}//namespace Adelus

#endif
