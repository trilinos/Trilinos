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

#ifndef __ADELUS_FACTOR_HPP__
#define __ADELUS_FACTOR_HPP__

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "Adelus_defines.h"
#include "Adelus_macros.h"
#include "Adelus_pcomm.hpp"
#include "Adelus_mytime.hpp"

#include "Kokkos_Core.hpp"
#include "KokkosBlas1_scal.hpp"
#include "BlasWrapper_copy.hpp"
#include "KokkosBlas1_iamax.hpp"
#include "KokkosBlas3_gemm.hpp"

extern int myrow;
extern int mycol;
extern int me;	               // processor id information
extern int nprocs_row;         // num of procs to which a row is assigned
extern int nprocs_col;         // num of procs to which a col is assigned
extern int nrows_matrix;       // number of rows in the matrix
extern int ncols_matrix;       // number of cols in the matrix
extern int my_rows;            // num of rows I own
extern int my_cols;            // num of cols I own
extern int my_rhs;             // num of right hand side I own
extern int blksz;              // block size for BLAS 3 operations

#define LUSTATUSINT 64

extern MPI_Comm col_comm;

//  Message tags
#define LUPIVOTTYPE (1<<13)
#define LUCOLTYPE (1<<14)
#define LUROWTYPE (1<<15)
#define LUPIVROWTYPE (1<<18)
#define LUSENDTYPE (1<<19)

#define rowplus(I) proc_num(mesh_row(me),(mesh_col(me)+(I) < nprocs_row) ? mesh_col(me)+(I) : mesh_col(me)+(I)-nprocs_row)
#define rowminus(I) rowplus(nprocs_row - (I))
#define MAXDIST 1

namespace Adelus {

template<class ZDView, class ViewType1D, class ViewType2D, class ViewIntType1D>
inline
void factor(ZDView& ZV,                    // matrix and rhs
            ViewType2D& col1_view,         // col used for updating a col
            ViewType2D& row1_view,         // diagonal row
            ViewType1D& row2_view,         // pivot row
            ViewType1D& row3_view,         // temporary vector for rows
            ViewIntType1D& pivot_vec_view) // vector storing list of pivot rows
{
  typedef typename ZDView::value_type value_type;
#ifdef PRINT_STATUS
  typedef typename ZDView::device_type::execution_space execution_space;
  typedef typename ZDView::device_type::memory_space memory_space;
#endif
#if defined(CUDA_HOST_PINNED_MPI) && defined(KOKKOS_ENABLE_CUDA)
  typedef Kokkos::View<value_type*, Kokkos::LayoutLeft, Kokkos::CudaHostPinnedSpace> View1DHostPinnType;//CudaHostPinnedSpace
#endif
  
  int j,k;               // loop counters

  struct pivot_type {    // a stucture for storing pivot info
    value_type entry;    //   pivot entry
    value_type current;  //   current row entry
    int row;             //   pivot row number
  } pivot;               // pivot info and some temporary pivot info

  constexpr int one = 1; 
  value_type d_one     = static_cast<value_type>(1.0);
  value_type d_min_one = static_cast<value_type>(-1.0);
  value_type d_zero    = static_cast<value_type>(0.0);

  int colcnt,cols_in_blk_owned; // number of columns stored for BLAS 3 ops
  int c_owner,r_owner,pivot_owner;
  int col_len,row_len,length,row_size,rows_used,cols_used;
  int rel_lpivot_row=0,lpivot_row;

  double pivot_mag;
  value_type invpiv;

  int gpivot_row; // make sure that this is well aligned

  int cur_col_i, cur_col_j, cur_row_i, cur_row_j, act_col_i, act_row_j, update_i, update_j;
  int sav_col_i, sav_col_j, sav_piv_row_i, sav_piv_row_j, act_piv_row_i, piv_row_i;
  int cur_col1_row_i, piv_col1_row_i;

  int ringdist,rdist;
  long type,bytes;

#if defined(DREAL) || defined(ZCPLX)
  struct {
    double val;
    int proc;
  } pivot_in,pivot_out;
#else
  struct {
    float val;
    int proc;
  } pivot_in,pivot_out;
#endif

  value_type entry,current;
  int row;
  int numprocs;

#ifdef GET_TIMING
  double updatetime,colupdtime,rowupdtime,scaltime;
  double xpivmsgtime,bcastpivstime,bcastpivrtime,bcastcolstime,bcastcolrtime,bcastrowtime,sendrowtime,recvrowtime;
  double copycoltime,copyrowtime,copyrow1time,copypivrowtime,copypivrow1time,pivotswaptime;
  double t1,t2;
  double msgtime,copytime,dgemmtime,totalfactortime;
  double iamaxtime,getlocalpivtime,localpivtime;
#if defined(CUDA_HOST_PINNED_MPI) && defined(KOKKOS_ENABLE_CUDA)
  double copyhostpinnedtime;
#endif
#endif

  MPI_Request msgrequest;
  MPI_Status msgstatus;

  // Distribution for the matrix on me

  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  if ( (numprocs/nprocs_row) * nprocs_row != numprocs ) {
     if (me == 0) {
       printf("nprocs_row must go into numprocs perfectly!\n");
       printf("Try a different value of nprocs_row.\n");
     }
     MPI_Barrier(MPI_COMM_WORLD);
     exit(0);
  }

#ifdef PRINT_STATUS
  printf("Rank %i -- factor() Begin LU factorization, execution_space %s, memory_space %s\n", me, typeid(execution_space).name(), typeid(memory_space).name());
#ifdef USE_DEEPCOPY
  printf("Rank %i -- factor() -- use Kokkos::deep_copy, KokkosBlas::iamax\n", me);
#else
  printf("Rank %i -- factor() -- use KokkosBlas::copy, KokkosBlas::iamax\n", me);
#endif
#endif
  
#ifdef GET_TIMING
  t2 = MPI_Wtime();
#endif
  colcnt = 0;                       // number of column's currently saved for update
  col_len = my_rows;                // length of column in remaining local matrix
 
  cur_col_i=0; cur_col_j=0;         // location of first column in remaining local matrix   
  sav_col_i=0; sav_col_j=0;         // location to store next active column   
  act_col_i=0;                      // location of matrix of columns being saved for gemm update
  
  row_len = my_cols + my_rhs;       // length of row in local matrix including rhs's
                                    
  rows_used = 0;                    // haven't used any local rows yet
  cols_used = 0;

  cur_row_i=0; cur_row_j=0;         // location of first row in local matrix
  cur_col1_row_i=0;                 // location of first row in col1 matrix
  act_row_j = 0;                    // location of matrix of rows being saved for gemm update
  sav_piv_row_i=0; sav_piv_row_j=0; // location for next row being saved for gemm update
  update_i=0; update_j=0;           // location of remaining local matrix

#ifdef GET_TIMING
  xpivmsgtime=bcastpivstime=bcastpivrtime=bcastcolstime=bcastcolrtime=bcastrowtime=sendrowtime=recvrowtime=0.0;
  copycoltime=copyrowtime=copyrow1time=copypivrowtime=copypivrow1time=pivotswaptime=0.0;
  updatetime=colupdtime=rowupdtime=scaltime=0.0;
  iamaxtime=getlocalpivtime=localpivtime=0.0;
#if defined(CUDA_HOST_PINNED_MPI) && defined(KOKKOS_ENABLE_CUDA)
  copyhostpinnedtime=0.0;
#endif
#endif

#ifdef PRINT_STATUS
  if (me == 0) {
    printf("Pivots sent separate from column \n");
    printf("Send to %4d neighbors \n",MAXDIST);
    printf("Attempt to do row work while waiting for column to arrive \n");
#ifdef OVERLAP
    printf("Column updates within block done one column at a time \n");
#endif
  }
#endif

#if defined(CUDA_HOST_PINNED_MPI) && defined(KOKKOS_ENABLE_CUDA)
  View1DHostPinnType h_coltmp( "h_coltmp", my_rows );
  View1DHostPinnType h_row2  ( "h_row2",   my_cols + blksz + nrhs );
  View1DHostPinnType h_row3  ( "h_row3",   my_cols + blksz + nrhs );
#endif

  //Kokkos::fence();
  for (j=0; j<ncols_matrix; j++) {
    c_owner = col_owner(j); r_owner = row_owner(j);
    ringdist = mesh_col(me) - mesh_col(c_owner);
    if (ringdist < 0) ringdist += nprocs_row;
    if (me == c_owner) {
      if (col_len > 0) {

#ifndef OVERLAP
        if (colcnt != 0){ // update current column with saved columns
#ifdef GET_TIMING
          t1 = MPI_Wtime();
#endif
          auto act_col_view = subview(col1_view,Kokkos::make_pair(act_col_i, act_col_i+col_len),Kokkos::make_pair(0, colcnt));
          auto act_row_view = subview(row1_view,Kokkos::make_pair(0, colcnt),                   Kokkos::make_pair(act_row_j, act_row_j+one));
          auto cur_col_view = subview(ZV,       Kokkos::make_pair(cur_col_i, cur_col_i+col_len),Kokkos::make_pair(cur_col_j, cur_col_j+one));
          KokkosBlas::gemm("N","N",d_min_one,
                           act_col_view,
                           act_row_view,
                           d_one,
                           cur_col_view);
#ifndef GET_TIMING
          if (numprocs > 1)
            Kokkos::fence();//Note: add Kokkos::fence() to guarantee synchronization if GET_TIMING not defined
#else
          Kokkos::fence();
          colupdtime += (MPI_Wtime()-t1);
#endif
        }
#endif

        // find maximum local pivot
#ifdef GET_TIMING
        t1 = MPI_Wtime();
#endif
        auto cur_col_view_1d = subview(ZV,Kokkos::make_pair(cur_col_i, cur_col_i+col_len),cur_col_j);
        rel_lpivot_row = KokkosBlas::iamax(cur_col_view_1d) - 1;
#ifdef GET_TIMING
        iamaxtime += (MPI_Wtime()-t1);
#endif
#ifdef GET_TIMING
        t1 = MPI_Wtime();
#endif

        Kokkos::View<value_type, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > pivot_entry(&pivot.entry);
        Kokkos::deep_copy(pivot_entry, subview(cur_col_view_1d,rel_lpivot_row));

        pivot.row = lrow_to_grow(rows_used+rel_lpivot_row);
        if (me == r_owner) {
          Kokkos::View<value_type, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > pivot_current(&pivot.current);
          Kokkos::deep_copy(pivot_current, subview(cur_col_view_1d,0));
        }
        else
          pivot.current = d_zero;
#ifdef GET_TIMING
        getlocalpivtime += (MPI_Wtime()-t1);
#endif
      } else {
#ifdef GET_TIMING
        t1 = MPI_Wtime();
#endif
        pivot.row = 0;  pivot.entry = d_zero;  pivot.current = d_zero;
#ifdef GET_TIMING
        iamaxtime += 0.0;
        getlocalpivtime += (MPI_Wtime()-t1);
#endif
      }
#ifdef GET_TIMING
      t1 = MPI_Wtime();
#endif

      // Pivot for column j on me

      entry   = pivot.entry;
      row     = pivot.row;
      current = pivot.current;
      pivot_in.val  = abs(pivot.entry);
      pivot_in.proc = me;
      bytes = sizeof(ADELUS_DATA_TYPE);
#ifdef DEBUG
      printf("Node %d: pivot val %g, pivot row %d \n",me,pivot_in.val,pivot.row);
#endif

      // Exchange to find global pivot value
      MPI_Allreduce(&pivot_in,&pivot_out,1,ADELUS_MPI_DATA_TYPE2,MPI_MAXLOC,col_comm);
      // Changed the mesh_row argument
      MPI_Bcast(&current,bytes,MPI_BYTE,mesh_row(r_owner),col_comm);
      MPI_Barrier(col_comm);
      pivot.current = current;
      MPI_Bcast(&row,1,MPI_INT,mesh_row(pivot_out.proc),col_comm);
      MPI_Barrier(col_comm);
      pivot.row = row;
      MPI_Bcast(&entry,bytes,MPI_BYTE,mesh_row(pivot_out.proc),col_comm);
      MPI_Barrier(col_comm);
      pivot.entry = entry;

#ifdef GET_TIMING
      xpivmsgtime += (MPI_Wtime()-t1);
#endif

      gpivot_row = pivot.row;
      pivot_mag = abs(pivot.entry);
      if (pivot_mag == 0.0) {
        printf("Node %d error -- zero pivot found in column %d -- exiting\n",me,j);
        return; 
      }

      // divide everything including the diagonal by the pivot entry

#ifdef GET_TIMING
      t1 = MPI_Wtime();
#endif
      invpiv = d_one/pivot.entry;
      auto cur_col_view_1d = subview(ZV, Kokkos::make_pair(cur_col_i, cur_col_i+col_len),cur_col_j);
      KokkosBlas::scal(cur_col_view_1d,invpiv,cur_col_view_1d);
#ifndef GET_TIMING
      if (numprocs > 1)
        Kokkos::fence();//Note: add Kokkos::fence() to guarantee synchronization if GET_TIMING not defined
#else
      Kokkos::fence();
      scaltime += (MPI_Wtime()-t1);
#endif

#ifdef GET_TIMING
      t1 = MPI_Wtime();
#endif
      // restore the pivot entry
      // swap pivot and current row in current column
      if (me == row_owner(gpivot_row)) {
        value_type pivot_swap_scalar = pivot.current * invpiv;
        Kokkos::View<value_type, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > pivot_swap(&pivot_swap_scalar);
        Kokkos::deep_copy(subview(ZV,cur_col_i+rel_lpivot_row,cur_col_j), pivot_swap);
      }
      if (me == r_owner) {
        Kokkos::View<value_type, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > pivot_entry(&pivot.entry);
        Kokkos::deep_copy(subview(cur_col_view_1d,0), pivot_entry);
      }
#ifdef GET_TIMING
      pivotswaptime += (MPI_Wtime()-t1);
#endif

#ifdef GET_TIMING
      t1 = MPI_Wtime();
#endif
      auto sav_col_view_1d = subview(col1_view,Kokkos::make_pair(sav_col_i, sav_col_i+col_len),sav_col_j);
#ifdef USE_DEEPCOPY
      Kokkos::deep_copy (sav_col_view_1d, cur_col_view_1d);//copy: cur_col_view_1d --> sav_col_view_1d
#else
      BlasWrapper::copy (cur_col_view_1d, sav_col_view_1d);//copy: cur_col_view_1d --> sav_col_view_1d
#endif
#ifndef GET_TIMING
      if (numprocs > 1)
        Kokkos::fence();//Note: add Kokkos::fence() to guarantee synchronization if GET_TIMING not defined
#else
      Kokkos::fence();
      copycoltime += (MPI_Wtime()-t1);
#endif

#ifdef GET_TIMING
      t1 = MPI_Wtime();
#endif
      // send column and pivot down one's row for column j
      for (rdist = 1;rdist <= MAXDIST;rdist++){
        if (rowplus(rdist) == c_owner) break;
        bytes = sizeof(gpivot_row);
        MPI_Send(&gpivot_row,bytes,MPI_BYTE,rowplus(rdist),LUPIVROWTYPE+j,MPI_COMM_WORLD);
      }
#ifdef GET_TIMING
      bcastpivstime += (MPI_Wtime()-t1);
#endif

#if defined(CUDA_HOST_PINNED_MPI) && defined(KOKKOS_ENABLE_CUDA)
#ifdef GET_TIMING
      t1 = MPI_Wtime();
#endif
      Kokkos::deep_copy(subview(h_coltmp, Kokkos::make_pair(0, col_len)), subview(col1_view, Kokkos::make_pair(sav_col_i, sav_col_i+col_len), sav_col_j));
#ifdef GET_TIMING
      copyhostpinnedtime += (MPI_Wtime()-t1);
#endif
#endif

#ifdef GET_TIMING
      t1 = MPI_Wtime();
#endif
      for (rdist = 1;rdist <= MAXDIST;rdist++) {
        if (rowplus(rdist) == c_owner) break;
        bytes=sizeof(ADELUS_DATA_TYPE)*col_len;
#if defined(CUDA_HOST_PINNED_MPI) && defined(KOKKOS_ENABLE_CUDA)
        MPI_Send(h_coltmp.data(),bytes,MPI_BYTE,rowplus(rdist),LUROWTYPE+j,MPI_COMM_WORLD);
#else //CUDA-aware MPI
        MPI_Send(col1_view.data()+sav_col_j*col1_view.stride(1)+sav_col_i,bytes,MPI_BYTE,rowplus(rdist),LUROWTYPE+j,MPI_COMM_WORLD);
#endif
      }
#ifdef GET_TIMING
      bcastcolstime += (MPI_Wtime()-t1);
#endif

      // if own active column don't include it in row work anymore

      row_len--; 
      update_j++;
      cur_col_j++;
      cur_row_j++;
      act_row_j++;
      sav_piv_row_j++;
      cols_used++;
    }
    else {

      // recv column and pivot

      bytes=col_len*sizeof(ADELUS_DATA_TYPE);
#if defined(CUDA_HOST_PINNED_MPI) && defined(KOKKOS_ENABLE_CUDA)
      MPI_Irecv(h_coltmp.data(),bytes,MPI_BYTE,MPI_ANY_SOURCE,LUROWTYPE+j,MPI_COMM_WORLD,&msgrequest);
#else //CUDA-aware MPI
      MPI_Irecv(col1_view.data()+sav_col_j*col1_view.stride(1)+sav_col_i,bytes,MPI_BYTE,
                MPI_ANY_SOURCE,LUROWTYPE+j,MPI_COMM_WORLD,&msgrequest);
#endif

#ifdef GET_TIMING
      t1 = MPI_Wtime();
#endif
      bytes = 0; type = LUPIVROWTYPE+j;
      bytes=4;
      bytes = sizeof(gpivot_row);
      MPI_Recv(&gpivot_row,bytes,MPI_BYTE,MPI_ANY_SOURCE,type,MPI_COMM_WORLD,&msgstatus);
#ifdef GET_TIMING
      bcastpivrtime += (MPI_Wtime()-t1);
#endif

      // if necessary forward column and pivot

      if ((ringdist % MAXDIST) == 0) {
#ifdef GET_TIMING
        t1 = MPI_Wtime();
#endif
        for (rdist = 1;rdist <= MAXDIST;rdist++) {
          if (rowplus(rdist) == c_owner) break;
          bytes = sizeof(gpivot_row);
          MPI_Send(&gpivot_row,bytes,MPI_BYTE,rowplus(rdist),LUPIVROWTYPE+j,MPI_COMM_WORLD);
        }
#ifdef GET_TIMING
        bcastpivstime += (MPI_Wtime()-t1);
#endif

#ifdef GET_TIMING
        t1 = MPI_Wtime();
#endif
        MPI_Wait(&msgrequest,&msgstatus);
#ifdef GET_TIMING
        bcastcolrtime += (MPI_Wtime()-t1);
#endif

#if defined(CUDA_HOST_PINNED_MPI) && defined(KOKKOS_ENABLE_CUDA)
#ifdef GET_TIMING
        t1 = MPI_Wtime();
#endif
        bytes=col_len*sizeof(ADELUS_DATA_TYPE);
        Kokkos::deep_copy(subview(col1_view, Kokkos::make_pair(sav_col_i, sav_col_i+col_len), sav_col_j), subview(h_coltmp, Kokkos::make_pair(0, col_len)));
#ifdef GET_TIMING
        copyhostpinnedtime += (MPI_Wtime()-t1);
#endif
#endif

#ifdef GET_TIMING
        t1 = MPI_Wtime();
#endif
        for (rdist = 1;rdist <= MAXDIST;rdist++) {
          if (rowplus(rdist) == c_owner) break;
          bytes=col_len*sizeof(ADELUS_DATA_TYPE);
#if defined(CUDA_HOST_PINNED_MPI) && defined(KOKKOS_ENABLE_CUDA)
          MPI_Send(h_coltmp.data(),bytes,MPI_BYTE,rowplus(rdist),LUROWTYPE+j,MPI_COMM_WORLD);
#else //CUDA-aware MPI
          MPI_Send(col1_view.data()+sav_col_j*col1_view.stride(1)+sav_col_i,bytes,MPI_BYTE,rowplus(rdist),LUROWTYPE+j,MPI_COMM_WORLD);
#endif
        }
#ifdef GET_TIMING
        bcastcolstime += (MPI_Wtime()-t1);
#endif
      }
    }
    pivot_owner = row_owner(gpivot_row); lpivot_row = grow_to_lrow(gpivot_row);
    act_piv_row_i = lpivot_row;
    piv_row_i = cur_col_i + (lpivot_row - rows_used);
    piv_col1_row_i = act_col_i + (lpivot_row - rows_used);
    row_size = (row_len + colcnt)*sizeof(ADELUS_DATA_TYPE);

    // send current row to owner of pivot row, skip this if pivot row is current row

    if (gpivot_row != j) {
      if (me == r_owner) { // own current row so pack it up

        // first take current row then portion in stored active columns
        // stored portion is at act_cur_row_ptr with col1's stride

#ifdef GET_TIMING
        t1 = MPI_Wtime();
#endif
        // copy from matrix and col1

        auto cur_row_view_1d = subview(ZV, cur_row_i, Kokkos::make_pair(cur_row_j, cur_row_j+row_len));
        auto sav_row_view_1d = subview(row2_view, Kokkos::make_pair(0, 0+row_len));  //note: sav_row_ptr = row2;
#ifdef USE_DEEPCOPY
        Kokkos::deep_copy (sav_row_view_1d, cur_row_view_1d);//copy: cur_row_view_1d --> sav_row_view_1d 
#else
        BlasWrapper::copy (cur_row_view_1d, sav_row_view_1d);//copy: cur_row_view_1d --> sav_row_view_1d
#endif
        auto cur_col1_row_view_1d = subview(col1_view, cur_col1_row_i, Kokkos::make_pair(0, 0+colcnt));
        auto sav_row_plus_view_1d = subview(row2_view, Kokkos::make_pair(row_len, row_len+colcnt));  //note: sav_row_ptr = row2;
#ifdef USE_DEEPCOPY
        Kokkos::deep_copy (sav_row_plus_view_1d, cur_col1_row_view_1d);//copy: cur_col1_row_view_1d --> sav_row_plus_view_1d 
#else
        BlasWrapper::copy (cur_col1_row_view_1d, sav_row_plus_view_1d);//copy: cur_col1_row_view_1d --> sav_row_plus_view_1d
#endif
#ifndef GET_TIMING
        if (numprocs > 1)
          Kokkos::fence();//Note: add Kokkos::fence() to guarantee synchronization if GET_TIMING not defined
#else
        Kokkos::fence();
        copyrowtime += (MPI_Wtime()-t1);
#endif

      }
    }

    // update pivot row and save in active row matrix

    if (me == pivot_owner) {
      if (colcnt > 0) {
#ifdef GET_TIMING
        t1 = MPI_Wtime();
#endif
#ifdef OVERLAP
        // don't need to update columns which were already updated

        cols_in_blk_owned = 0;
        for (k = 1; k < (blksz - colcnt); k++)
          if (j+k < ncols_matrix)
            if (me == col_owner(j+k)) cols_in_blk_owned++;
        length = row_len - cols_in_blk_owned;
        if (length > 0) {
          auto act_piv_row_view = subview(col1_view,Kokkos::make_pair(act_piv_row_i, act_piv_row_i+one),Kokkos::make_pair(0, 0+colcnt));
          auto act_row_view     = subview(row1_view,Kokkos::make_pair(0, colcnt), Kokkos::make_pair(act_row_j+cols_in_blk_owned, act_row_j+cols_in_blk_owned+length));
          auto piv_row_view     = subview(ZV, Kokkos::make_pair(piv_row_i, piv_row_i+one),Kokkos::make_pair(cur_col_j+cols_in_blk_owned, cur_col_j+cols_in_blk_owned+length));
          KokkosBlas::gemm("N","N",d_min_one,
                           act_piv_row_view,
                           act_row_view,
                           d_one,
                           piv_row_view);
#ifndef GET_TIMING
          if (numprocs > 1)
            Kokkos::fence();//Note: add Kokkos::fence() to guarantee synchronization if GET_TIMING not defined
#endif
        }
#else
        auto act_piv_row_view = subview(col1_view,Kokkos::make_pair(act_piv_row_i, act_piv_row_i+one),Kokkos::make_pair(0, 0+colcnt));
        auto act_row_view     = subview(row1_view,Kokkos::make_pair(0, colcnt), Kokkos::make_pair(act_row_j+cols_in_blk_owned, act_row_j+cols_in_blk_owned+length));
        auto piv_row_view     = subview(ZV, Kokkos::make_pair(piv_row_i, piv_row_i+one),Kokkos::make_pair(cur_col_j+cols_in_blk_owned, cur_col_j+cols_in_blk_owned+length));
        KokkosBlas::gemm("N","N",d_min_one,
                         act_piv_row_view,
                         act_row_view,
                         d_one,
                         piv_row_view);
#ifndef GET_TIMING
        if (numprocs > 1)
          Kokkos::fence();//Note: add Kokkos::fence() to guarantee synchronization if GET_TIMING not defined
#endif
#endif
#ifdef GET_TIMING
        Kokkos::fence();        
        rowupdtime += (MPI_Wtime()-t1);
#endif
      }

      // copy pivot row to temp holder

#ifdef GET_TIMING
      t1 = MPI_Wtime();
#endif

      // copy from matrix and col1

      auto piv_row_view_1d  = subview(ZV, piv_row_i, Kokkos::make_pair(cur_col_j, cur_col_j+row_len));
      auto temp_row_view_1d = subview(row3_view, Kokkos::make_pair(0, 0+row_len));  //note: temp_row_ptr -> row3.view;
#ifdef USE_DEEPCOPY
      Kokkos::deep_copy (temp_row_view_1d, piv_row_view_1d);//copy: piv_row_view_1d --> temp_row_view_1d
#else
      BlasWrapper::copy (piv_row_view_1d, temp_row_view_1d);//copy: piv_row_view_1d --> temp_row_view_1d
#endif

      auto piv_col1_row_view_1d  = subview(col1_view, piv_col1_row_i, Kokkos::make_pair(0, 0+colcnt));
      auto temp_row_plus_view_1d = subview(row3_view, Kokkos::make_pair(row_len, row_len+colcnt));  //note: temp_row_ptr -> row3.view;
#ifdef USE_DEEPCOPY
      Kokkos::deep_copy (temp_row_plus_view_1d, piv_col1_row_view_1d);//copy: piv_col1_row_view_1d --> temp_row_plus_view_1d 
#else
      BlasWrapper::copy (piv_col1_row_view_1d, temp_row_plus_view_1d);//copy: piv_col1_row_view_1d --> temp_row_plus_view_1d
#endif
#ifndef GET_TIMING
      if (numprocs > 1)
        Kokkos::fence();//Note: add Kokkos::fence() to guarantee synchronization if GET_TIMING not defined
#else
      Kokkos::fence();
      copypivrowtime += (MPI_Wtime()-t1);
#endif
    }

    // broadcast pivot row

#if defined(CUDA_HOST_PINNED_MPI) && defined(KOKKOS_ENABLE_CUDA)
#ifdef GET_TIMING
    t1 = MPI_Wtime();
#endif
    Kokkos::deep_copy(subview(h_row3, Kokkos::make_pair(0, row_len + colcnt)), subview(row3_view, Kokkos::make_pair(0, row_len + colcnt)));
#ifdef GET_TIMING
    copyhostpinnedtime += (MPI_Wtime()-t1);
#endif
#endif

#ifdef GET_TIMING
    t1 = MPI_Wtime();
#endif
    bytes=sizeof(ADELUS_DATA_TYPE)*row_size ;
#if defined(CUDA_HOST_PINNED_MPI) && defined(KOKKOS_ENABLE_CUDA)
    MPI_Bcast(reinterpret_cast<char *>(h_row3.data()), row_size, MPI_CHAR, mesh_row(pivot_owner), col_comm);
    MPI_Barrier(col_comm);
#else //CUDA-aware MPI -- Note: Looks like MPI_Bcast is still working well with device (cuda) pointers (and faster than using cuda host pinned memory) for 2 nodes
    MPI_Bcast(reinterpret_cast<char *>(row3_view.data()), row_size, MPI_CHAR, mesh_row(pivot_owner), col_comm);
    MPI_Barrier(col_comm);
#endif
#ifdef GET_TIMING
    bcastrowtime += (MPI_Wtime()-t1);
#endif

#if defined(CUDA_HOST_PINNED_MPI) && defined(KOKKOS_ENABLE_CUDA)
#ifdef GET_TIMING
    t1 = MPI_Wtime();
#endif
    Kokkos::deep_copy(subview(row3_view, Kokkos::make_pair(0, row_len + colcnt)), subview(h_row3, Kokkos::make_pair(0, row_len + colcnt)));
#ifdef GET_TIMING
    copyhostpinnedtime += (MPI_Wtime()-t1);
#endif
#endif

#ifdef GET_TIMING
    t1 = MPI_Wtime();
#endif

    auto temp_row_view_1d    = subview(row3_view, Kokkos::make_pair(0, 0+row_len));  //note: temp_row_ptr -> row3.view;
    auto sav_piv_row_view_1d = subview(row1_view, sav_piv_row_i, Kokkos::make_pair(sav_piv_row_j, sav_piv_row_j+row_len));
#ifdef USE_DEEPCOPY
    Kokkos::deep_copy (sav_piv_row_view_1d, temp_row_view_1d);//copy: temp_row_view_1d --> sav_piv_row_view_1d
#else
    BlasWrapper::copy (temp_row_view_1d, sav_piv_row_view_1d);//copy: temp_row_view_1d --> sav_piv_row_view_1d
#endif
#ifndef GET_TIMING
    if (numprocs > 1)
      Kokkos::fence();//Note: add Kokkos::fence() to guarantee synchronization if GET_TIMING not defined
#else
    Kokkos::fence();
    copypivrowtime += (MPI_Wtime()-t1);
#endif

    if (gpivot_row != j) {
      if (me != pivot_owner && me == r_owner) {
#if defined(CUDA_HOST_PINNED_MPI) && defined(KOKKOS_ENABLE_CUDA)
#ifdef GET_TIMING
        t1 = MPI_Wtime();
#endif
        Kokkos::deep_copy(subview(h_row2, Kokkos::make_pair(0, row_len + colcnt)), subview(row2_view, Kokkos::make_pair(0, row_len + colcnt)));
#ifdef GET_TIMING
        copyhostpinnedtime += (MPI_Wtime()-t1);
#endif
#endif

#ifdef GET_TIMING
        t1 = MPI_Wtime();
#endif
        bytes=(row_len+colcnt)*sizeof(ADELUS_DATA_TYPE);
#if defined(CUDA_HOST_PINNED_MPI) && defined(KOKKOS_ENABLE_CUDA)
        MPI_Send(h_row2.data(),bytes,MPI_BYTE,pivot_owner,LUSENDTYPE+j,MPI_COMM_WORLD);
#else //CUDA-aware MPI
        MPI_Send(row2_view.data(),bytes,MPI_BYTE,pivot_owner,LUSENDTYPE+j,MPI_COMM_WORLD);
#endif
#ifdef GET_TIMING
        sendrowtime += (MPI_Wtime()-t1);
#endif
      }
      if (me == pivot_owner) {
        // receive top row and copy into pivot row

#ifdef GET_TIMING
        t1 = MPI_Wtime();
#endif
        if (me != r_owner) {
          bytes=(row_len+colcnt)*sizeof(ADELUS_DATA_TYPE);
#if defined(CUDA_HOST_PINNED_MPI) && defined(KOKKOS_ENABLE_CUDA)
          MPI_Recv(h_row2.data(),bytes,MPI_BYTE,r_owner,LUSENDTYPE+j,MPI_COMM_WORLD,&msgstatus);
#else //CUDA-aware MPI
          MPI_Recv(row2_view.data(),bytes,MPI_BYTE,r_owner,LUSENDTYPE+j,MPI_COMM_WORLD,&msgstatus);
#endif
        }
#ifdef GET_TIMING
        recvrowtime += (MPI_Wtime()-t1);
#endif

        if (me != r_owner) {
#if defined(CUDA_HOST_PINNED_MPI) && defined(KOKKOS_ENABLE_CUDA)
#ifdef GET_TIMING
          t1 = MPI_Wtime();
#endif
          Kokkos::deep_copy(subview(row2_view, Kokkos::make_pair(0, row_len + colcnt)), subview(h_row2, Kokkos::make_pair(0, row_len + colcnt)));
#ifdef GET_TIMING
          copyhostpinnedtime += (MPI_Wtime()-t1);
#endif
#endif
        }

#ifdef GET_TIMING
        t1 = MPI_Wtime();
#endif
        // copy from matrix and col1

        auto sav_row_view_1d = subview(row2_view, Kokkos::make_pair(0, 0+row_len));  //note: sav_row_ptr = row2;
        auto piv_row_view_1d = subview(ZV, piv_row_i, Kokkos::make_pair(cur_col_j, cur_col_j+row_len));
#ifdef USE_DEEPCOPY
        Kokkos::deep_copy (piv_row_view_1d, sav_row_view_1d);//copy: sav_row_view_1d --> piv_row_view_1d
#else
        BlasWrapper::copy (sav_row_view_1d, piv_row_view_1d);//copy: sav_row_view_1d --> piv_row_view_1d
#endif

        auto sav_row_plus_view_1d = subview(row2_view, Kokkos::make_pair(row_len, row_len+colcnt));  //note: sav_row_ptr = row2;
        auto piv_col1_row_view_1d  = subview(col1_view, piv_col1_row_i, Kokkos::make_pair(0, 0+colcnt));
#ifdef USE_DEEPCOPY
        Kokkos::deep_copy (piv_col1_row_view_1d, sav_row_plus_view_1d);//copy: sav_row_plus_view_1d --> piv_col1_row_view_1d
#else
        BlasWrapper::copy (sav_row_plus_view_1d, piv_col1_row_view_1d);//copy: sav_row_plus_view_1d --> piv_col1_row_view_1d
#endif
#ifndef GET_TIMING
        if (numprocs > 1)
          Kokkos::fence();//Note: add Kokkos::fence() to guarantee synchronization if GET_TIMING not defined
#else
        Kokkos::fence();
        copyrow1time += (MPI_Wtime()-t1);
#endif
      }
      if (me == r_owner) { /* copy pivot row into current row */
#ifdef GET_TIMING
        t1 = MPI_Wtime();
#endif

        // copy from matrix and col1

        auto tmp_row_view_1d = subview(row3_view, Kokkos::make_pair(0, 0+row_len));  //note: temp_row_ptr -> row3.view;
        auto cur_row_view_1d = subview(ZV, cur_row_i, Kokkos::make_pair(cur_row_j, cur_row_j+row_len));
#ifdef USE_DEEPCOPY
        Kokkos::deep_copy (cur_row_view_1d, tmp_row_view_1d);//copy: tmp_row_view_1d --> cur_row_view_1d
#else
        BlasWrapper::copy (tmp_row_view_1d, cur_row_view_1d);//copy: tmp_row_view_1d --> cur_row_view_1d
#endif
        auto temp_row_plus_view_1d = subview(row3_view, Kokkos::make_pair(row_len, row_len+colcnt));  //note: temp_row_ptr -> row3.view;
        auto cur_col1_row_view_1d  = subview(col1_view, cur_col1_row_i, Kokkos::make_pair(0, 0+colcnt));
#ifdef USE_DEEPCOPY
        Kokkos::deep_copy (cur_col1_row_view_1d, temp_row_plus_view_1d);//copy: temp_row_plus_view_1d --> cur_col1_row_view_1d
#else
        BlasWrapper::copy (temp_row_plus_view_1d, cur_col1_row_view_1d);//copy: temp_row_plus_view_1d --> cur_col1_row_view_1d
#endif
#ifndef GET_TIMING
        if (numprocs > 1)
          Kokkos::fence();//Note: add Kokkos::fence() to guarantee synchronization if GET_TIMING not defined
#else
        Kokkos::fence();
        copypivrow1time += (MPI_Wtime()-t1);
#endif
      }
    }
    if ((me != c_owner) && ((ringdist % MAXDIST) != 0)) {
#ifdef GET_TIMING
      t1 = MPI_Wtime();
#endif
      MPI_Wait(&msgrequest,&msgstatus);
#ifdef GET_TIMING
      bcastcolrtime += (MPI_Wtime()-t1);
#endif

#if defined(CUDA_HOST_PINNED_MPI) && defined(KOKKOS_ENABLE_CUDA)
#ifdef GET_TIMING
      t1 = MPI_Wtime();
#endif
      Kokkos::deep_copy(subview(col1_view, Kokkos::make_pair(sav_col_i, sav_col_i+col_len), sav_col_j), subview(h_coltmp, Kokkos::make_pair(0, col_len)));
#ifdef GET_TIMING
      copyhostpinnedtime += (MPI_Wtime()-t1);
#endif
#endif
    }
    // saved this active row and column so get ready for next ones

    if (me == r_owner) { // finished with this row so update all column pointers
      col_len--; rows_used++;
      update_i++;
      cur_row_i++;
      cur_col1_row_i++;
      cur_col_i++;
      sav_col_i++;
      act_col_i++;
    }
    colcnt++;

#ifdef OVERLAP
    cols_in_blk_owned = 0;
    for (k = 1; k <= (blksz - colcnt); k++)
      if (j+k < ncols_matrix)
        if (me == col_owner(j+k)) cols_in_blk_owned++;
    if (cols_in_blk_owned > 0){ /* update current column with latest column */
#ifdef GET_TIMING
      t1 = MPI_Wtime();
#endif
      auto sav_col_view     = subview(col1_view,Kokkos::make_pair(sav_col_i, sav_col_i+col_len),Kokkos::make_pair(sav_col_j, sav_col_j+one));
      auto sav_piv_row_view = subview(row1_view,Kokkos::make_pair(sav_piv_row_i, sav_piv_row_i+one), Kokkos::make_pair(sav_piv_row_j, sav_piv_row_j+cols_in_blk_owned));
      auto cur_col_view     = subview(ZV, Kokkos::make_pair(cur_col_i, cur_col_i+col_len),Kokkos::make_pair(cur_col_j, cur_col_j+cols_in_blk_owned));
      KokkosBlas::gemm("N","N",d_min_one,
                       sav_col_view,
                       sav_piv_row_view,
                       d_one,
                       cur_col_view);
#ifndef GET_TIMING
      if (numprocs > 1)
        Kokkos::fence();//Note: add Kokkos::fence() to guarantee synchronization if GET_TIMING not defined
#else
      Kokkos::fence();
      colupdtime += (MPI_Wtime()-t1);
#endif
    }
#endif
 
    sav_col_j++; 
    sav_piv_row_i++;

    // if we have saved up enough columns, we do the outer product update

    if (colcnt == blksz)
      if (j != ncols_matrix-1){
#ifdef GET_TIMING
        t1 = MPI_Wtime();
#endif
        auto act_col_view=subview(col1_view,Kokkos::make_pair(act_col_i, act_col_i+col_len),Kokkos::make_pair(0, colcnt));
        auto act_row_view=subview(row1_view,Kokkos::make_pair(0, colcnt),                   Kokkos::make_pair(act_row_j, act_row_j+row_len));
        auto update_view =subview(ZV,       Kokkos::make_pair(update_i, update_i+col_len),  Kokkos::make_pair(update_j, update_j+row_len));
        KokkosBlas::gemm("N","N",d_min_one,
                         act_col_view,
                         act_row_view,
                         d_one,
                         update_view);
#ifndef GET_TIMING
        if (numprocs > 1)
          Kokkos::fence();//Note: add Kokkos::fence() to guarantee synchronization if GET_TIMING not defined
#else
        Kokkos::fence();
        updatetime += (MPI_Wtime()-t1);
#endif
        // reset active matrix pointers

        colcnt = 0; 
        act_col_i=rows_used; sav_col_i=rows_used; sav_col_j=0;
        act_row_j = 0; sav_piv_row_i=sav_piv_row_j=0;
     }
//#ifdef PRINT_STATUS
//    if (((j%1000) == 0) && (me == 0)) {
//      fprintf(stderr," Column %d completed\n",j);
//    }
//#endif
  }

#ifdef GET_TIMING
  totalfactortime = MPI_Wtime() - t2;

  localpivtime = iamaxtime+getlocalpivtime;
  msgtime      = xpivmsgtime+bcastpivstime+bcastpivrtime+bcastcolstime+bcastcolrtime+bcastrowtime+sendrowtime+recvrowtime;
  copytime     = pivotswaptime+copycoltime+copyrowtime+copyrow1time+copypivrowtime+copypivrow1time;
  dgemmtime    = updatetime+colupdtime+rowupdtime+scaltime;
#ifdef ADELUS_SHOW_TIMING_DETAILS
  showtime("Time to do iamax",&iamaxtime);
  showtime("Time to get local pivot",&getlocalpivtime);
#endif
  showtime("Total finding local pivot time",&localpivtime);
  double tmp = 100*localpivtime/totalfactortime;
  showtime("Percent finding local pivot time",&tmp);
#ifdef ADELUS_SHOW_TIMING_DETAILS
  showtime("Time to xchgpivot",&xpivmsgtime);
  showtime("Time to do send in bcast pivot",&bcastpivstime);
  showtime("Time to do recv in bcast pivot",&bcastpivrtime);
  tmp = bcastpivrtime+bcastpivstime;
  showtime("Time to do bcast pivot",&tmp);
  showtime("Time to do send in bcast cur col",&bcastcolstime);
  showtime("Time to do recv bcast cur col",&bcastcolrtime);
  tmp = bcastcolrtime+bcastcolstime;
  showtime("Time to do bcast cur col",&tmp);
  tmp = bcastcolrtime+bcastcolstime+bcastpivrtime+bcastpivstime;
  showtime("Time to do bcast cur col and pivot",&tmp);
  showtime("Time to bcast piv row",&bcastrowtime);
  showtime("Time to send cur row",&sendrowtime);
  showtime("Time to recv cur row",&recvrowtime);
#endif
  showtime("Total msg passing time",&msgtime);
  tmp = 100*msgtime/totalfactortime;
  showtime("Percent msg passing time",&tmp);
#if defined(CUDA_HOST_PINNED_MPI) && defined(KOKKOS_ENABLE_CUDA)
  showtime("Total copy between host pinned mem and dev mem time",&copyhostpinnedtime); 
  tmp = 100*copyhostpinnedtime/totalfactortime;
  showtime("Percent copy between host pinned mem and dev mem time",&tmp);  
#endif
#ifdef ADELUS_SHOW_TIMING_DETAILS
  showtime("Time to swap pivot",&pivotswaptime);
  showtime("Time to copy cur col",&copycoltime);
  showtime("Time to copy cur row to sav row",&copyrowtime);
  showtime("Time to copy piv row to sav piv",&copypivrowtime);
  showtime("Time to copy sav row to cur row",&copyrow1time);
  showtime("Time to copy sav piv  to piv row",&copypivrow1time);
#endif
  showtime("Total copying time",&copytime);
  tmp = 100*copytime/totalfactortime;
  showtime("Percent copying time",&tmp);
#ifdef ADELUS_SHOW_TIMING_DETAILS
  showtime("Time to scale cur col",&scaltime);
  showtime("Time to update cur col",&colupdtime);
  showtime("Time to update piv row",&rowupdtime);
#endif
  showtime("Time to update matrix",&updatetime);
  showtime("Total update time",&dgemmtime);
  tmp = 100*dgemmtime/totalfactortime;
  showtime("Percent update time",&tmp);
  showtime("Total time in factor",&totalfactortime);
#endif
}

}//namespace Adelus

#endif
