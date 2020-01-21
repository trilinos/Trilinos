/*
//@HEADER
// ************************************************************************
//
//               Pliris: Parallel Dense Solver Package
//                 Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
#ifndef __SOLVE_HPP__
#define __SOLVE_HPP__

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "defines.h"
//#include "BLAS_prototypes_cpp.hpp"
#include "macros.h"
#include "pcomm.h"
#include "mytime.hpp"

#include "Kokkos_Core.hpp"
#include "KokkosBlas3_gemm.hpp"

extern int me;

extern int ncols_matrix;       /* number of cols in the matrix */

extern int nprocs_col;		/* num of procs to which a col is assigned */
extern int nprocs_row;		/* num of procs to which a row is assigned */

extern int my_first_col;        /* proc position in a col */
extern int my_first_row;	/* proc position in a row */

extern int my_rows;		/* num of rows I own */
extern int my_cols;            /* num of cols I own */

extern int nrhs;                /* number of right hand sides */
extern int my_rhs;              /* number of right hand sides that I own */

//extern int rhs_blksz;

extern MPI_Comm col_comm;


#define SOSTATUSINT 32768

/*  Previous message tags  
#define SOCOLTYPE (1<<26)
#define SOROWTYPE (1<<27)
#define SOHSTYPE (SOCOLTYPE + SOROWTYPE)
#define PERMTYPE ((1 << 25) + (1 << 26))

*/

/* New message tags  */

#define SOCOLTYPE (1<<14)
#define SOROWTYPE (1<<15)
#define SOHSTYPE (SOCOLTYPE + SOROWTYPE)
//#define PERMTYPE ((1 << 14) + (1 << 15))

/*  Customized elimination on the rhs that I own */	
//diag_mag = ABS_VAL(*(ptr2+end_row-1));
//for (j1=0; j1<n_rhs_this; j1++) {
//  ptr3 = rhs + j1*my_rows + end_row - 1;
//  ptr4 = row1 + j1;
//  DIVIDE(*ptr3,*(ptr2+end_row-1),diag_mag,*ptr4);
//  *ptr3 = *ptr4;
//}	  
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
void back_solve6(ZDView& ZV)
{
  typedef typename ZDView::value_type value_type;
  typedef typename ZDView::device_type::execution_space execution_space;
  typedef typename ZDView::device_type::memory_space memory_space;
  typedef Kokkos::View<value_type**, Kokkos::LayoutLeft, memory_space> ViewMatrixType;

#if defined(CUDA_HOST_PINNED_MPI) && defined(KOKKOS_ENABLE_CUDA)
  typedef Kokkos::View<value_type**, Kokkos::LayoutLeft, Kokkos::CudaHostPinnedSpace> View2DHostPinnType;//CudaHostPinnedSpace
#endif

  int  j;                      /* loop counters */

  int end_row;                  /* row num to end column operations */

  //DATA_TYPE *ptr2;             /* running ptrs into mat and update cols */
  //DATA_TYPE *ptr3, *ptr4;       /* running ptrs to diagonal and pivot rows */
  //DATA_TYPE *ptr5;              /* running ptr into matrix */

  int bytes[16];                /* number of bytes in messages */
  int root;                     /* root processor for fanout */
  int type[16];                 /* mesage type for messages */
  int dest[16];                 /* dest for message sends */

  int one = 1;                  /* constant for the daxpy routines */

  //DATA_TYPE d_one = CONST_ONE;  /* constant for the daxpy routines */
  //DATA_TYPE d_min_one = CONST_MINUS_ONE; /* constant for the daxpy routines */
  value_type d_one = 1.0;/* constant for the daxpy routines */
  value_type d_min_one = -1.0;/* constant for the daxpy routines */
  
  int /*j1,*/ j2;

  //double diag_mag;              /* magnitude of matrix diagonal */

  int n_rhs_this;               /* num rhs that I currently own */

  //char transA = 'N';            /* all dgemm's don't transpose matrix A  */
  //char transB = 'N';            /* nor transpose matrix B */

  int col_offset;               /* which processor starts the pipeline */
  int my_pos;                   /* my position in the new linup */
  int extra;                    /* extra loop to realign data after pipeline */
  int act_col;                  /* act this column (that I own) */
  int on_col;                   /* on this collection of rhs's */
  int global_col;               /* global col number for act_col */
  int max_bytes;                /* max number of bytes of rhs I can receive */

  int my_col_id, my_row_id, id_temp;
  int dest_right, dest_left;

  int blas_length;

#ifdef GET_TIMING
  double t1,t2;
  double allocviewtime,eliminaterhstime,bcastrowtime,updrhstime,sendrhstime,recvrhstime,copyrhstime;
  double totalsolvetime;
#if defined(CUDA_HOST_PINNED_MPI) && defined(KOKKOS_ENABLE_CUDA)
  double copyhostpinnedtime;
#endif
#endif

  MPI_Request msgrequest;
  MPI_Status msgstatus;

  //DATA_TYPE *fseg;
  //DATA_TYPE *rhs;
  //DATA_TYPE *row1;
  //fseg =  reinterpret_cast<DATA_TYPE *>(ZV.data());
  //rhs  =  reinterpret_cast<DATA_TYPE *>(ZV.data()+my_rows*my_cols);
  //row1 =  reinterpret_cast<DATA_TYPE *>(row1_view.data());

#ifdef GET_TIMING
  t2 = MPI_Wtime();
#endif

  /* set the block size for backsolve */

  my_col_id = mesh_col(me);
  my_row_id = mesh_row(me);

  id_temp = my_col_id + 1;
  if (id_temp >= nprocs_row) id_temp = 0;
  dest_right = proc_num(my_row_id,id_temp);

  id_temp = my_col_id - 1;
  if (id_temp < 0) id_temp = nprocs_row-1;
  dest_left = proc_num(my_row_id,id_temp);

  /* set j2 to be first column in last group of columns */
  //rhs_blksz=1;
  max_bytes = nrhs/nprocs_row;                     //printf("Rank %i -- back_solve6() max_bytes (1st) %d\n",me, max_bytes);
  if (nrhs%nprocs_row > 0) max_bytes++;            //printf("Rank %i -- back_solve6() max_bytes (2nd) %d\n",me, max_bytes);
  max_bytes = max_bytes*sizeof(DATA_TYPE)*my_rows; //printf("Rank %i -- back_solve6() max_bytes (3rd) %d\n",me, max_bytes);

#ifdef GET_TIMING
  allocviewtime=eliminaterhstime=bcastrowtime=updrhstime=sendrhstime=recvrhstime=copyrhstime=0.0;
#if defined(CUDA_HOST_PINNED_MPI) && defined(KOKKOS_ENABLE_CUDA)
  copyhostpinnedtime=0.0;
#endif
#endif

  //DATA_TYPE *row2 = (DATA_TYPE *) malloc( max_bytes);
  //if (row2 == NULL) {
  // fprintf(stderr, "Node %d: Out of memory\n", me);
  // exit(-1);
  //}

#ifdef GET_TIMING
  t1 = MPI_Wtime();
#endif
  ViewMatrixType row1( "row1", one, nrhs );   /* row1: diagonal row (temp variables) */
  ViewMatrixType row2( "row2", my_rows, max_bytes/sizeof(DATA_TYPE)/my_rows );

#if defined(CUDA_HOST_PINNED_MPI) && defined(KOKKOS_ENABLE_CUDA)
  View2DHostPinnType h_row1( "h_row1", one, nrhs );
  View2DHostPinnType h_row2( "h_row2", my_rows, max_bytes/sizeof(DATA_TYPE)/my_rows );
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

  printf("Rank %i -- back_solve6() Begin back solve, execution_space %s, memory_space %s\n",me, typeid(execution_space).name(), typeid(memory_space).name());
  
  for (j = j2; j >= 1-nprocs_row-extra; j--) {

    if ((j+nprocs_row-1 >= 0) && (n_rhs_this > 0)) {

      if ((act_col < my_cols) && (act_col >= 0)) {

        global_col = act_col*nprocs_row + my_first_col;

        end_row = global_col/nprocs_col;
        if (my_first_row <= global_col%nprocs_col) ++end_row;

        //ptr5 = fseg + act_col*my_rows;
                 
        /* do an elimination step on the rhs that I own */

        //ptr2 = ptr5;
        auto ptr2_view = subview(ZV, end_row-1, Kokkos::ALL());//*(ptr2+end_row-1)

        root = row_owner(global_col);

        if (me == root) {
#ifdef GET_TIMING
          t1 = MPI_Wtime();
#endif
          //diag_mag = ABS_VAL(*(ptr2+end_row-1));
          //for (j1=0; j1<n_rhs_this; j1++) {
          //  ptr3 = rhs + j1*my_rows + end_row - 1;
          //  ptr4 = reinterpret_cast<DATA_TYPE *>(row1.data()) + j1;
          //  DIVIDE(*ptr3,*(ptr2+end_row-1),diag_mag,*ptr4);
          //  *ptr3 = *ptr4;
          //}
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
        bytes[0] = n_rhs_this*sizeof(DATA_TYPE);
        type[0]  = SOCOLTYPE+j;

        //MPI_Bcast((char *) row1, bytes[0], MPI_CHAR, mesh_row(root), col_comm);
//#if defined(CUDA_HOST_PINNED_MPI) && defined(KOKKOS_ENABLE_CUDA)
//        MPI_Bcast(reinterpret_cast<char *>(h_row1.data()), bytes[0], MPI_CHAR, mesh_row(root), col_comm);
//#else //CUDA-aware MPI -- Note: Looks like MPI_Bcast is still working well with device (cuda) pointers (and faster than using cuda host pinned memory)
        MPI_Bcast(reinterpret_cast<char *>(row1.data()), bytes[0], MPI_CHAR, mesh_row(root), col_comm);		
//#endif
        /* added this barrier for CPLANT operation  */

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
       /* Changed XGEMM_ to XGEMM   removed all &   */

        //XGEMM_(&transA, &transB, &end_row, &n_rhs_this, &one, &d_min_one,
        //       ptr5, &my_rows,
        //       row1, &one, &d_one,
        //       rhs, &my_rows);
        //XGEMM_(&transA, &transB, &end_row, &n_rhs_this, &one, &d_min_one,
        //       reinterpret_cast<DATA_TYPE *>(ZV.data()+my_rows*act_col), &my_rows,
        //       reinterpret_cast<DATA_TYPE *>(row1.data()), &one, &d_one,
        //       reinterpret_cast<DATA_TYPE *>(ZV.data()+my_rows*my_cols), &my_rows);

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

    if (j != 1-nprocs_row-extra) {
      dest[0] = dest_right;
      if (me != dest[0]) {
        bytes[0] = max_bytes;
        type[0]  = SOROWTYPE+j;

        //MPI_Irecv((char *) row2, bytes[0], MPI_CHAR, MPI_ANY_SOURCE,type[0],MPI_COMM_WORLD,&msgrequest);
#if defined(CUDA_HOST_PINNED_MPI) && defined(KOKKOS_ENABLE_CUDA)
        MPI_Irecv(reinterpret_cast<char *>(h_row2.data()), bytes[0], MPI_CHAR, MPI_ANY_SOURCE, type[0], MPI_COMM_WORLD, &msgrequest);
#else //CUDA-aware MPI
        MPI_Irecv(reinterpret_cast<char *>(  row2.data()), bytes[0], MPI_CHAR, MPI_ANY_SOURCE, type[0], MPI_COMM_WORLD, &msgrequest);
#endif

        n_rhs_this = bytes[0]/sizeof(DATA_TYPE)/my_rows;

#if defined(CUDA_HOST_PINNED_MPI) && defined(KOKKOS_ENABLE_CUDA)
#ifdef GET_TIMING
        t1 = MPI_Wtime();
#endif
        Kokkos::deep_copy(subview(h_rhs, Kokkos::ALL(), Kokkos::make_pair(0, n_rhs_this)), subview(ZV, Kokkos::ALL(), Kokkos::make_pair(my_cols, my_cols+n_rhs_this)));
#ifdef GET_TIMING
        copyhostpinnedtime += (MPI_Wtime()-t1);
#endif
#endif

#ifdef GET_TIMING
        t1 = MPI_Wtime();
#endif
        dest[1]  = dest_left;
        bytes[1] = n_rhs_this * sizeof(DATA_TYPE) * my_rows;
        type[1]  = SOROWTYPE+j;

        //MPI_Send((char *) rhs, bytes[1], MPI_CHAR, dest[1], type[1],MPI_COMM_WORLD);
#if defined(CUDA_HOST_PINNED_MPI) && defined(KOKKOS_ENABLE_CUDA)
        MPI_Send(reinterpret_cast<char *>(h_rhs.data()), bytes[1], MPI_CHAR, dest[1], type[1], MPI_COMM_WORLD);
#else //CUDA-aware MPI
        MPI_Send(reinterpret_cast<char *>(ZV.data()+my_rows*my_cols), bytes[1], MPI_CHAR, dest[1], type[1], MPI_COMM_WORLD);
#endif
#ifdef GET_TIMING
        sendrhstime += (MPI_Wtime()-t1);
#endif

#ifdef GET_TIMING
        t1 = MPI_Wtime();
#endif
        MPI_Wait(&msgrequest,&msgstatus);
#ifdef GET_TIMING
        recvrhstime += (MPI_Wtime()-t1);
#endif

#if defined(CUDA_HOST_PINNED_MPI) && defined(KOKKOS_ENABLE_CUDA)
#ifdef GET_TIMING
        t1 = MPI_Wtime();
#endif
        Kokkos::deep_copy(row2,h_row2);
#ifdef GET_TIMING
        copyhostpinnedtime += (MPI_Wtime()-t1);
#endif
#endif

#ifdef GET_TIMING
        t1 = MPI_Wtime();
#endif
        // Copy row2 -> rhs
        blas_length = n_rhs_this*my_rows;
        //XCOPY(blas_length,row2,one,rhs,one);
        //XCOPY(blas_length,reinterpret_cast<DATA_TYPE *>(row2.data()),one,reinterpret_cast<DATA_TYPE *>(ZV.data()+my_rows*my_cols),one);
        //deep_copy(subview(ZV,   Kokkos::ALL(), Kokkos::make_pair(my_cols, my_cols+n_rhs_this)),
        //          subview(row2, Kokkos::ALL(), Kokkos::make_pair(0, n_rhs_this)));//deep_copy is slower than BLAS XCOPY
#ifdef KOKKOS_ENABLE_CUDA//Use memcpy for now
        cudaMemcpy(reinterpret_cast<DATA_TYPE *>(ZV.data()+my_rows*my_cols), reinterpret_cast<DATA_TYPE *>(row2.data()), blas_length*sizeof(DATA_TYPE), cudaMemcpyDeviceToDevice);
#else
        memcpy(reinterpret_cast<DATA_TYPE *>(ZV.data()+my_rows*my_cols), reinterpret_cast<DATA_TYPE *>(row2.data()), blas_length*sizeof(DATA_TYPE));
#endif
#ifdef GET_TIMING
        copyrhstime += (MPI_Wtime()-t1);
#endif
      }
      on_col++;
      if (on_col >= nprocs_row) {
        on_col = 0;
        act_col--;
      }
    }
  }
  /* free(row2);  */
#ifdef GET_TIMING
  totalsolvetime = MPI_Wtime() - t2;
#endif
#ifdef GET_TIMING
  showtime("Time to alloc view",&allocviewtime);
  showtime("Time to eliminate rhs",&eliminaterhstime);
  showtime("Time to bcast temp row",&bcastrowtime);
  showtime("Time to update rhs",&updrhstime);
  showtime("Time to send in rhs",&sendrhstime);
  showtime("Time to recv rhs",&recvrhstime);
#if defined(CUDA_HOST_PINNED_MPI) && defined(KOKKOS_ENABLE_CUDA)
  showtime("Time to copy host pinned mem <--> dev mem",&copyhostpinnedtime);   
#endif
  showtime("Time to copy rhs",&copyrhstime);
  showtime("Total time in solve",&totalsolvetime);
#endif
}

//void collect_vector(DATA_TYPE *vec)
//{
//  int j, k;
//
//  int start_col;
//  int end_row;
//
//  int bytes;
//  int dest;
//  int type;
//
//  MPI_Status msgstatus;
//
//  for (j=0; j<ncols_matrix; j++) {
//    if (me == col_owner(j)) {
//      start_col = (j) / nprocs_row;
//      if (my_first_col < (j)%nprocs_row) ++start_col;
//
//      if (j%rhs_blksz == 0) {
//        for (k=0; k<rhs_blksz; k++) {
//          end_row = (j+k)/nprocs_col;
//          if (my_first_row <= (j+k)%nprocs_col) ++end_row;
//          if (j+k < ncols_matrix) {
//            if (me == row_owner(j+k)) {
//              dest = col_owner(0);
//              if (me != dest) {
//                bytes = sizeof(DATA_TYPE);
//                type = PERMTYPE + j + k;
//                MPI_Send((vec+end_row-1),bytes,MPI_BYTE,
//                   dest,type,MPI_COMM_WORLD);
//              }
//            }
//          }
//        }
//      }
//    }
//    if (me == col_owner(0)) {
//      if (me == row_owner(j)) {
//        end_row = (j)/nprocs_col;
//        if (my_first_row <= (j)%nprocs_col) ++end_row;
//        dest = col_owner((j/rhs_blksz)*rhs_blksz);
//        if (me != dest) {
//          bytes = sizeof(DATA_TYPE);
//          type = PERMTYPE + j;
//          MPI_Recv((vec+end_row-1),bytes,MPI_BYTE,dest,
//                   type,MPI_COMM_WORLD,&msgstatus);
//        }
//      }
//    }
//  }
//}
#endif
