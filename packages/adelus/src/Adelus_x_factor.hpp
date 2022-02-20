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

#ifndef __ADELUS_XLU_HPP__
#define __ADELUS_XLU_HPP__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Adelus_defines.h"
#include "mpi.h"
#include "Adelus_vars.hpp"
#include "Adelus_macros.h"
#include "Adelus_block.h"
#include "Adelus_factor.hpp"
#include "Adelus_perm_mat.hpp"
#include "Adelus_pcomm.hpp"
#include "Adelus_mytime.hpp"
#include "Kokkos_Core.hpp"

#ifdef ADELUS_HAVE_TIME_MONITOR
#include "Teuchos_TimeMonitor.hpp"
#endif

namespace Adelus {

template<class ZDView, class IDView>
inline
void lu_(ZDView& ZV, int *matrix_size, int *num_procsr, IDView& permute, double *secs)
{
#ifdef ADELUS_HAVE_TIME_MONITOR
  using Teuchos::TimeMonitor;
#endif

  using value_type      = typename ZDView::value_type;
#ifdef PRINT_STATUS
  using execution_space = typename ZDView::device_type::execution_space;
#endif
  using memory_space    = typename ZDView::device_type::memory_space;

  double run_secs;              // time (in secs) during which the prog ran
  double tsecs;                 // intermediate storage of timing info
  int totmem;

  // Determine who I am (me ) and the total number of nodes (nprocs_cube)
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs_cube);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  nrows_matrix = *matrix_size;
  ncols_matrix = *matrix_size;
  nprocs_row   = *num_procsr;

  totmem=0;                      // Initialize the total memory used
  nprocs_col = nprocs_cube/nprocs_row;
  max_procs = (nprocs_row < nprocs_col) ? nprocs_col : nprocs_row;

  // Set up communicators for rows and columns
  myrow = mesh_row(me);
  mycol = mesh_col(me);

  MPI_Comm_split(MPI_COMM_WORLD,myrow,mycol,&row_comm);

  MPI_Comm_split(MPI_COMM_WORLD,mycol,myrow,&col_comm);

  // Distribution for the matrix on me
  my_first_col = mesh_col(me);
  my_first_row = mesh_row(me);

  my_rows = nrows_matrix / nprocs_col;
  if (my_first_row < nrows_matrix % nprocs_col)
    ++my_rows;
  my_cols = ncols_matrix / nprocs_row;
  if (my_first_col < ncols_matrix % nprocs_row)
    ++my_cols;

  // blksz parameter must be set
  blksz = DEFBLKSZ;

#ifdef PRINT_STATUS
  printf("Rank %i -- factor_() Begin LU with blksz %d, value_type %s, execution_space %s, memory_space %s\n", me, blksz, typeid(value_type).name(), typeid(execution_space).name(), typeid(memory_space).name());
#endif

  // Allocate arrays for factor/solve
  typedef Kokkos::View<value_type*,  Kokkos::LayoutLeft, memory_space> ViewType1D;
  typedef Kokkos::View<value_type**, Kokkos::LayoutLeft, memory_space> ViewType2D;

  totmem += (blksz) * (my_rows) * sizeof(ADELUS_DATA_TYPE);             //col1_view
  totmem += blksz * (my_cols + blksz + nrhs) * sizeof(ADELUS_DATA_TYPE);//row1_view
  totmem += (my_cols + blksz + nrhs) * sizeof(ADELUS_DATA_TYPE);        //row2_view
  totmem += (my_cols + blksz + nrhs) * sizeof(ADELUS_DATA_TYPE);        //row3_view
  totmem += my_cols * sizeof(int);                                      //lpiv_view
  
  ViewType2D  col1_view ( "col1_view", my_rows, blksz );
  ViewType2D  row1_view ( "row1_view", blksz, my_cols + blksz + nrhs );
  ViewType1D  row2_view ( "row2_view", my_cols + blksz + nrhs );
  ViewType1D  row3_view ( "row3_view", my_cols + blksz + nrhs );
  IDView      lpiv_view ( "lpiv_view", my_cols );

  {
  // Factor the system

  tsecs = get_seconds(0.0);

  initcomm();

#ifdef PRINT_STATUS
  printf("OpenMP or Cuda: Rank %i -- factor() starts ...\n", me);
#endif
#ifdef ADELUS_HAVE_TIME_MONITOR
  {
    TimeMonitor t(*TimeMonitor::getNewTimer("Adelus: factor"));
#endif
    factor(ZV,
           col1_view,
           row1_view,
           row2_view, 
           row3_view, 
           lpiv_view);
#ifdef ADELUS_HAVE_TIME_MONITOR
  }
#endif

  // Permute the lower triangular matrix
  //NOTE: Currently doing matrix permutation in host memory
  typename ZDView::HostMirror h_ZV = Kokkos::create_mirror_view( ZV );
  Kokkos::deep_copy (h_ZV, ZV);

#ifdef ADELUS_HAVE_TIME_MONITOR
  {
    TimeMonitor t(*TimeMonitor::getNewTimer("Adelus: matrix permutation"));
#endif
    permute_mat(h_ZV, lpiv_view, permute);
#ifdef ADELUS_HAVE_TIME_MONITOR
  }
#endif

  Kokkos::deep_copy (ZV, h_ZV);

  tsecs = get_seconds(tsecs);

  run_secs = (double) tsecs;
  
  *secs = run_secs;
  showtime("Total time in Factor (inl. matrix permutation)",&run_secs);
  }
}

}//namespace Adelus

#endif
