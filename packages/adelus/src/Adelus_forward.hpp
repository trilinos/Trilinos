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

#ifndef __ADELUS_FORWARD_HPP__
#define __ADELUS_FORWARD_HPP__

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "Adelus_defines.h"
#include "Adelus_macros.h"
#include "Adelus_pcomm.hpp"
#include "Adelus_mytime.hpp"
#include "Kokkos_Core.hpp"

namespace Adelus {

template<class ZView, class RHSView>
inline
void forward(ZView& Z, RHSView& RHS)
{
  //NOTE: Currently assume that Z and RHS reside in host memory, and 
  //                            there is only a single RHS vector
  using value_type = typename ZView::value_type ;
#ifdef PRINT_STATUS
  using execution_space = typename ZView::device_type::execution_space ;
  using memory_space    = typename ZView::device_type::memory_space ;
#endif
  using ViewVectorType  =  Kokkos::View<value_type*, Kokkos::LayoutLeft, Kokkos::HostSpace>;

  int rhs_col;     // torus-wrap column containing the rhs
  int k_row;       // torus-wrap row corresponding to kth global row
  int k_col;       // torus-wrap column corresponding to kth global col
  int istart;      // Starting row index for pivot column
  int count_row;   // dummy index

  value_type ck;   // rhs corresponding to current column of the backsubstitution
  ViewVectorType piv_col( "piv_col", my_rows ); // portion of pivot column I am sending

  MPI_Request msgrequest;
  MPI_Status msgstatus;

#ifdef PRINT_STATUS
  printf("Rank %i -- forward() Begin forward solve with myrow %d, mycol %d, nprocs_row %d, nprocs_col %d, nrows_matrix %d, ncols_matrix %d, my_rows %d, my_cols %d, my_rhs %d, nrhs %d, value_type %s, execution_space %s, memory_space %s\n", me, myrow, mycol, nprocs_row, nprocs_col, nrows_matrix, ncols_matrix, my_rows, my_cols, my_rhs, nrhs, typeid(value_type).name(), typeid(execution_space).name(), typeid(memory_space).name());
#endif

#ifdef GET_TIMING
  double t1, fwdsolvetime;
  t1 = MPI_Wtime();
#endif

  // Perform the Forward Substitution:
  rhs_col = 0;
  for (int k=0; k<= nrows_matrix-2; k++) {
    k_row=k%nprocs_col;
    k_col=k%nprocs_row;
    istart = (k+1-myrow)/nprocs_col;
    if (istart * nprocs_col < k+1-myrow) istart++;
    count_row = 0;
    for (int i=istart;i<=my_rows-1;i++) {
      piv_col(count_row)=Z(i,k/nprocs_row);
      count_row++;
    }
    if (mycol == rhs_col && myrow == k_row) ck = RHS(k/nprocs_col,0);
    if (mycol == rhs_col) {
      MPI_Irecv(reinterpret_cast<char *>(piv_col.data()), count_row*sizeof(ADELUS_DATA_TYPE), MPI_CHAR, k_col, 0, row_comm, &msgrequest);
    }
    if (mycol == k_col) {
      MPI_Send(reinterpret_cast<char *>(piv_col.data()), count_row*sizeof(ADELUS_DATA_TYPE), MPI_CHAR, rhs_col, 0, row_comm);
    }
    if (mycol == rhs_col) {
      MPI_Wait(&msgrequest,&msgstatus);
    }
    if (mycol == rhs_col) {
      MPI_Bcast((char *)(&ck),sizeof(ADELUS_DATA_TYPE),MPI_CHAR,k_row,col_comm);
      count_row=0;

      for (int i=istart;i<=my_rows-1;i++) {
        RHS(i,0) = RHS(i,0) - piv_col(count_row) * ck;
        count_row++;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }// end of for (k=0; k<= nrows_matrix-2; k++)

#ifdef GET_TIMING
  fwdsolvetime = MPI_Wtime() - t1;
  showtime("Total time in forward solve",&fwdsolvetime);
#endif
}

}//namespace Adelus

#endif
