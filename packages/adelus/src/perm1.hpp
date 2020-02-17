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

#ifndef __PERM1_HPP__
#define __PERM1_HPP__

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "defines.h"
#include "macros.h"
#include "mytime.hpp"

#include "Kokkos_Core.hpp"

#define PERMTYPE ((1 << 5) + (1 << 4))

extern int me;	               /* processor id information */
extern int nprocs_row;         /* num of procs to which a row is assigned */
extern int nprocs_col;         /* num of procs to which a col is assigned */
extern int nrows_matrix;       /* number of rows in the matrix */
extern int ncols_matrix;       /* number of cols in the matrix */
extern int my_rows;            /* num of rows I own */
extern int my_cols;            /* num of cols I own */
extern int myrow;
extern int mycol;
extern MPI_Comm col_comm;

namespace Adelus {

/*  Customized zcopy  */
template<class XView, class YView>
void zcopy_wr_local_index(int N, XView& X, YView& Y, int lidx) {
  Kokkos::parallel_for(Kokkos::RangePolicy<typename XView::device_type::execution_space>(0,N+1), KOKKOS_LAMBDA (const int i) {
    if (i<N)
      Y(i) = X(i);
    else
      Y(i) = (double)lidx + 0.1;
  });
}

/*  Customized zcopy  */
template<class XView, class YView>
void zcopy_ld_local_index(int N, XView& X, YView& Y) {
  Kokkos::parallel_for(Kokkos::RangePolicy<typename XView::device_type::execution_space>(0,N), KOKKOS_LAMBDA (const int i) {
#ifdef COMPLEX
    int lidx = (int)(X(N).real());
#else
    int lidx = (int)(X(N));
#endif
    Y(lidx,i) = X(i);
  });
}

/*  Permutes -- unwraps the torus-wrap for the solution
    using the communication buffer             */
template<class ZDView>
inline
void perm1_(ZDView& ZV, int *num_my_rhs) {

  int i;
  int my_rhs_;


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
  int change_recv;
  int next, next_r;
  int inc;

#ifdef GET_TIMING
  double t2;
  double totalpermtime;
#endif

  MPI_Request *msgrequest;

  int ptr1_idx;

#ifdef GET_TIMING
  t2 = MPI_Wtime();
#endif

  my_rhs_=*num_my_rhs;
  
  typedef typename ZDView::device_type::execution_space execution_space;
  typedef typename ZDView::device_type::memory_space memory_space;
  typedef Kokkos::View<Kokkos::complex<double>*, Kokkos::LayoutLeft, memory_space>  ViewVectorType;
					 
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
    change_nosend = 0;
    next = 0;
    change_send = 0;
    change_recv = 0;
    next_r = 0;

#ifdef PRINT_STATUS
    printf("Rank %i -- perm1_() Begin permutation, execution_space %s, memory_space %s\n",me,typeid(execution_space).name(),typeid(memory_space).name());
#endif

    int col_comm_rank, col_comm_size;

    MPI_Comm_rank(col_comm, &col_comm_rank);
    MPI_Comm_size(col_comm, &col_comm_size);

    msgrequest = (MPI_Request *)malloc(col_comm_size * sizeof(MPI_Request));

    int j, dest_tmp, global_index_tmp, local_index_tmp, my_rows_max, recv_count;

    MPI_Allreduce( &my_rows, &my_rows_max, 1, MPI_INT, MPI_MAX, col_comm);

    for (i=0; i<my_rows_max; i++) {

      //Find out how many ranks will send data to me in this ith iteration
      //and what dest me will send data to 
      recv_count = 0;

      for (j=0; j<col_comm_size; j++) {//loop through ranks in the column communicator

        if (j == col_comm_rank) {//me
          if (i < my_rows) {//only do with my_rows
            global_index = my_first_row + i*nprocs_col;
            /* break down global index using torus wrap in row direction */
            dest = global_index%nprocs_row;
            local_index = global_index/nprocs_row;
            
            /* rebuild global index using block in row direction */
            if (dest < nprocs_row1) {
              global_index = dest*ncols_proc1 + local_index;
            } else {
              global_index = row_offset + (dest-nprocs_row1)*ncols_proc2
                             + local_index;
            }
            
            /* break down global index using blocks in the column direction */
            if (global_index < col_offset) {
              dest = global_index/nrows_proc1;
              local_index = global_index%nrows_proc1;
            } else {
              dest = (global_index - col_offset)/nrows_proc2 + nprocs_col1;
              local_index = (global_index - col_offset)%nrows_proc2;
            }
            
            dest = dest*nprocs_row + my_first_col;
          }
	    }
        else {//not me (me calculates for other ranks in the column communicator)
          int my_first_col_tmp = mesh_col(j*nprocs_row + mycol);//convert column communicator rank to global rank before applying macro
          int my_first_row_tmp = mesh_row(j*nprocs_row + mycol);//convert column communicator rank to global rank before applying macro
          int my_rows_tmp = nrows_matrix / nprocs_col;
          if (my_first_row_tmp < nrows_matrix % nprocs_col)
            my_rows_tmp++;
  
          if (i < my_rows_tmp) {//only do with the rows of this rank 
            global_index_tmp = my_first_row_tmp + i*nprocs_col;
            /* break down global index using torus wrap in row direction */
            dest_tmp = global_index_tmp%nprocs_row;
            local_index_tmp = global_index_tmp/nprocs_row;
            
            /* rebuild global index using block in row direction */
            if (dest_tmp < nprocs_row1) {
              global_index_tmp = dest_tmp*ncols_proc1 + local_index_tmp;
            } else {
              global_index_tmp = row_offset + (dest_tmp-nprocs_row1)*ncols_proc2
                             + local_index_tmp;
            }
            
            /* break down global index using blocks in the column direction */
            if (global_index_tmp < col_offset) {
              dest_tmp = global_index_tmp/nrows_proc1;
              local_index_tmp = global_index_tmp%nrows_proc1;
            } else {
              dest_tmp = (global_index_tmp - col_offset)/nrows_proc2 + nprocs_col1;
              local_index_tmp = (global_index_tmp - col_offset)%nrows_proc2;
            }
            
            dest_tmp = dest_tmp*nprocs_row + my_first_col_tmp;
            if (dest_tmp == me) recv_count++;
          }
        }

      }//End j loop

      bytes = (my_rhs_ + 1)*sizeof(DATA_TYPE);

      if (recv_count > 0) {
        //Receive data from other ranks
        for (j=0; j<recv_count; j++) {
            MPI_Irecv( (char *)(reinterpret_cast<DATA_TYPE *>(rhs_temp.data())+next_r), bytes, MPI_CHAR, 
                        MPI_ANY_SOURCE, MPI_ANY_TAG, col_comm, &msgrequest[j]);
            change_recv++;
            next_r = change_recv * (my_rhs_ + 1);
        }
      }

      if (i < my_rows) {//send data or change position within my_rows, otherwise only receive from other ranks  
        if ((local_index != i) || (dest != me)) {
        
          /* Check if I need to send the data or just change position */
        
          if( dest == me ) {
        
            auto sub_ZV          = subview(ZV,          ptr1_idx, Kokkos::ALL());
            auto sub_my_rhs_temp = subview(my_rhs_temp, Kokkos::make_pair(next, next+(my_rhs_ + 1)));     				
            zcopy_wr_local_index(my_rhs_, sub_ZV, sub_my_rhs_temp, local_index);

            change_nosend++;

            next = change_nosend * (my_rhs_ + 1);

          }
        
          if( dest !=me ) {
        
            auto sub_ZV = subview(ZV, ptr1_idx, Kokkos::ALL());     				
            zcopy_wr_local_index(my_rhs_, sub_ZV, temp_s, local_index);
        
            type = PERMTYPE + change_send;
        
            MPI_Send((char *)(reinterpret_cast<DATA_TYPE *>(temp_s.data())), bytes, MPI_CHAR, dest/nprocs_row,
                      type, col_comm);
            change_send++;

          }
        
        }
      }

      if (recv_count > 0) {
        MPI_Waitall( recv_count, msgrequest, MPI_STATUSES_IGNORE );
      }

      ptr1_idx++;

    }//End i loop

    /* Unpack changes from other processors  */

    next_r = 0;
    inc = 0;
    for (i = 0; i < change_recv; i++) {
      auto sub_rhs_temp = subview(rhs_temp, Kokkos::make_pair(next_r, next_r + (my_rhs_ + 1)));     				
      zcopy_ld_local_index(my_rhs_, sub_rhs_temp, ZV);

      inc++;
      next_r = inc * (my_rhs_+1);
    }

    /* Unpack my changes */
    next = 0;
    inc = 0;
    for (i = 0; i < change_nosend; i++) {
      auto sub_my_rhs_temp = subview(my_rhs_temp, Kokkos::make_pair(next, next+(my_rhs_ + 1)));     				
      zcopy_ld_local_index(my_rhs_, sub_my_rhs_temp, ZV);
      inc++;
      next = inc * (my_rhs_+1);
    }

    free(msgrequest);

  }

#ifdef GET_TIMING
  totalpermtime = MPI_Wtime() - t2;
#endif
#ifdef GET_TIMING
  showtime("Total time in perm",&totalpermtime);
#endif
}

}//namespace Adelus

#endif
