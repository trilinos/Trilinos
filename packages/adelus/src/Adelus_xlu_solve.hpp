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

#ifndef __ADELUS_XLUSOLVE_HPP__
#define __ADELUS_XLUSOLVE_HPP__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include "Kokkos_Core.hpp"
#include "Adelus_defines.h"
#include "Adelus_macros.h"
#include "Adelus_vars.hpp"
#include "Adelus_mytime.hpp"
#include "Adelus_solve.hpp"
#include "Adelus_factor.hpp"
#include "Adelus_perm1.hpp"

#ifdef ADELUS_HAVE_TIME_MONITOR
#include "Teuchos_TimeMonitor.hpp"
#endif

namespace Adelus {

template<class HandleType, class ZRHSViewType>
inline
void lusolve_(HandleType& ahandle, ZRHSViewType& ZRHS, double *secs)
{
#ifdef ADELUS_HAVE_TIME_MONITOR
  using Teuchos::TimeMonitor;
#endif

  using value_type      = typename ZRHSViewType::value_type;
#ifdef PRINT_STATUS
  using execution_space = typename ZRHSViewType::device_type::execution_space;
#endif
  using memory_space    = typename ZRHSViewType::device_type::memory_space;

  int blksz   = ahandle.get_blksz();
  int my_rows = ahandle.get_my_rows();
  int my_cols = ahandle.get_my_cols();
  int nrhs    = ahandle.get_nrhs();
  int my_rhs  = ahandle.get_my_rhs();

  double run_secs; // time (in secs) during which the prog ran
  double tsecs;    // intermediate storage of timing info
  int totmem = 0;  // Initialize the total memory used

#ifdef PRINT_STATUS
  printf("Rank %i -- lusolve_() Begin LU+Solve+Perm with blksz %d, value_type %s, execution_space %s, memory_space %s\n", ahandle.get_myrank(), blksz, typeid(value_type).name(), typeid(execution_space).name(), typeid(memory_space).name());
#endif

  // Allocate arrays for factor/solve
  typedef Kokkos::View<value_type*,  Kokkos::LayoutLeft, memory_space> ViewType1D;
  typedef Kokkos::View<value_type**, Kokkos::LayoutLeft, memory_space> ViewType2D;
  typedef Kokkos::View<int*, Kokkos::LayoutLeft, Kokkos::HostSpace> ViewIntType1D;

  totmem += (blksz) * (my_rows) * sizeof(ADELUS_DATA_TYPE);             //col1_view
  totmem += blksz * (my_cols + blksz + nrhs) * sizeof(ADELUS_DATA_TYPE);//row1_view
  totmem += (my_cols + blksz + nrhs) * sizeof(ADELUS_DATA_TYPE);        //row2_view
  totmem += (my_cols + blksz + nrhs) * sizeof(ADELUS_DATA_TYPE);        //row3_view
  totmem += my_cols * sizeof(int);                               //pivot_vec_view
  
  ViewType2D    col1_view      ( "col1_view",      my_rows, blksz );
  ViewType2D    row1_view      ( "row1_view",      blksz, my_cols + blksz + nrhs );
  ViewType1D    row2_view      ( "row2_view",      my_cols + blksz + nrhs );
  ViewType1D    row3_view      ( "row3_view",      my_cols + blksz + nrhs );  
  ViewIntType1D pivot_vec_view ( "pivot_vec_view", my_cols );

  {
  // Factor and Solve the system

  tsecs = get_seconds(0.0);

#ifdef PRINT_STATUS
  printf("OpenMP or Cuda: Rank %i -- factor() starts ...\n", ahandle.get_myrank());
#endif
#ifdef ADELUS_HAVE_TIME_MONITOR
  {
    TimeMonitor t(*TimeMonitor::getNewTimer("Adelus: factor"));
#endif
    factor(ahandle,
           ZRHS,
           col1_view,
           row1_view,
           row2_view, 
           row3_view, 
           pivot_vec_view,
           nrhs, my_rhs);
#ifdef ADELUS_HAVE_TIME_MONITOR
  }
#endif

  if (nrhs > 0) {
    auto Z   = subview(ZRHS, Kokkos::ALL(), Kokkos::make_pair(0, my_cols));
    auto RHS = subview(ZRHS, Kokkos::ALL(), Kokkos::make_pair(my_cols, my_cols + my_rhs + 6));

    // Perform the backsolve

#ifdef PRINT_STATUS
    printf("OpenMP or Cuda: Rank %i -- back_solve6() starts ...\n", ahandle.get_myrank());
#endif
#ifdef ADELUS_HAVE_TIME_MONITOR
    {
      TimeMonitor t(*TimeMonitor::getNewTimer("Adelus: backsolve"));
#endif
      back_solve6(ahandle, Z, RHS);
#ifdef ADELUS_HAVE_TIME_MONITOR
    }
#endif

    // Permute the results -- undo the torus map

#ifdef PRINT_STATUS
    printf("OpenMP or Cuda: Rank %i -- perm1_()(permute the results -- undo the torus map) starts ...\n", ahandle.get_myrank());
#endif
#ifdef ADELUS_HAVE_TIME_MONITOR
    {
      TimeMonitor t(*TimeMonitor::getNewTimer("Adelus: permutation"));
#endif
      perm1_(ahandle, RHS);
#ifdef ADELUS_HAVE_TIME_MONITOR
    }
#endif
  }

  tsecs = get_seconds(tsecs);

  run_secs = (double) tsecs;

  // Solve time secs

  *secs = run_secs;
  showtime(ahandle.get_comm_id(), ahandle.get_comm(), ahandle.get_myrank(), ahandle.get_nprocs_cube(),
           "Total time in Factor and Solve", &run_secs);
  }
}

}//namespace Adelus

#endif
