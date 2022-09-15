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
#include "mpi.h"
#include "Kokkos_Core.hpp"
#include "Adelus_defines.h"
#include "Adelus_macros.h"
#include "Adelus_vars.hpp"
#include "Adelus_mytime.hpp"
#include "Adelus_factor.hpp"
#include "Adelus_perm_mat.hpp"


#ifdef ADELUS_HAVE_TIME_MONITOR
#include "Teuchos_TimeMonitor.hpp"
#endif

namespace Adelus {

template<class HandleType, class ZViewType, class PViewType>
inline
void lu_(HandleType& ahandle, ZViewType& Z, PViewType& permute, double *secs)
{
#ifdef ADELUS_HAVE_TIME_MONITOR
  using Teuchos::TimeMonitor;
#endif

  using value_type      = typename ZViewType::value_type;
#ifdef PRINT_STATUS
  using execution_space = typename ZViewType::device_type::execution_space;
#endif
  using memory_space    = typename ZViewType::device_type::memory_space;

  int blksz   = ahandle.get_blksz();
  int my_rows = ahandle.get_my_rows();
  int my_cols = ahandle.get_my_cols();

  double run_secs; // time (in secs) during which the prog ran
  double tsecs;    // intermediate storage of timing info
  int totmem = 0;  // Initialize the total memory used
  
#ifdef PRINT_STATUS
  printf("Rank %i -- factor_() Begin LU with blksz %d, myrow %d, mycol %d, nprocs_row %d, nprocs_col %d, nrows_matrix %d, ncols_matrix %d, my_rows %d, my_cols %d, my_rhs %d, nrhs %d, value_type %s, execution_space %s, memory_space %s\n", ahandle.get_myrank(), blksz, ahandle.get_myrow(), ahandle.get_mycol(), ahandle.get_nprocs_row(), ahandle.get_nprocs_col(), ahandle.get_nrows_matrix(), ahandle.get_ncols_matrix(), my_rows, my_cols, ahandle.get_my_rhs(), ahandle.get_nrhs(), typeid(value_type).name(), typeid(execution_space).name(), typeid(memory_space).name());
#endif

  // Allocate arrays for factor
  using ViewType1D = Kokkos::View<value_type*,  Kokkos::LayoutLeft, memory_space>;
  using ViewType2D = Kokkos::View<value_type**, Kokkos::LayoutLeft, memory_space>;

  totmem += (blksz) * (my_rows) * sizeof(ADELUS_DATA_TYPE);          //col1_view
  totmem += blksz * (my_cols + blksz + 0) * sizeof(ADELUS_DATA_TYPE);//row1_view
  totmem += (my_cols + blksz + 0) * sizeof(ADELUS_DATA_TYPE);        //row2_view
  totmem += (my_cols + blksz + 0) * sizeof(ADELUS_DATA_TYPE);        //row3_view
  totmem += my_cols * sizeof(int);                                   //lpiv_view
  
  ViewType2D  col1_view ( "col1_view", my_rows, blksz );
  ViewType2D  row1_view ( "row1_view", blksz, my_cols + blksz + 0 );
  ViewType1D  row2_view ( "row2_view", my_cols + blksz + 0 );
  ViewType1D  row3_view ( "row3_view", my_cols + blksz + 0 );
  PViewType   lpiv_view ( "lpiv_view", my_cols );

  {
  // Factor the system

  tsecs = get_seconds(0.0);

#ifdef PRINT_STATUS
  printf("OpenMP or Cuda: Rank %i -- factor() starts ...\n", ahandle.get_myrank());
#endif
#ifdef ADELUS_HAVE_TIME_MONITOR
  {
    TimeMonitor t(*TimeMonitor::getNewTimer("Adelus: factor"));
#endif
    factor(ahandle,
           Z,
           col1_view,
           row1_view,
           row2_view, 
           row3_view, 
           lpiv_view,
           0, 0);
#ifdef ADELUS_HAVE_TIME_MONITOR
  }
#endif

  // Permute the lower triangular matrix
#ifdef ADELUS_HAVE_TIME_MONITOR
  {
    TimeMonitor t(*TimeMonitor::getNewTimer("Adelus: matrix permutation"));
#endif
#ifdef ADELUS_PERM_MAT_FORWARD_COPY_TO_HOST
    typename ZViewType::HostMirror h_Z = Kokkos::create_mirror_view( Z );
    Kokkos::deep_copy (h_Z, Z);
  
    permute_mat(ahandle, h_Z, lpiv_view, permute);

    Kokkos::deep_copy (Z, h_Z);
#else
    permute_mat(ahandle, Z, lpiv_view, permute);
#endif
#ifdef ADELUS_HAVE_TIME_MONITOR
  }
#endif

  tsecs = get_seconds(tsecs);

  run_secs = (double) tsecs;
  
  *secs = run_secs;
  showtime(ahandle.get_comm_id(), ahandle.get_comm(), ahandle.get_myrank(), ahandle.get_nprocs_cube(),
           "Total time in Factor (inl. matrix permutation)", &run_secs );
  }
}

}//namespace Adelus

#endif
