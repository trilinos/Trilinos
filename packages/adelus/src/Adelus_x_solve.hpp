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

#ifndef __ADELUS_XSOLVE_HPP__
#define __ADELUS_XSOLVE_HPP__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include "Kokkos_Core.hpp"
#include "Adelus_defines.h"
#include "Adelus_macros.h"
#include "Adelus_vars.hpp"
#include "Adelus_mytime.hpp"
#include "Adelus_perm_rhs.hpp"
#include "Adelus_forward.hpp"
#include "Adelus_solve.hpp"
#include "Adelus_perm1.hpp"

#ifdef ADELUS_HAVE_TIME_MONITOR
#include "Teuchos_TimeMonitor.hpp"
#endif

namespace Adelus {

template<class HandleType, class ZViewType, class RHSViewType, class PViewType>
inline
void solve_(HandleType& ahandle, ZViewType& Z, RHSViewType& RHS, PViewType& permute, double *secs)
{
#ifdef ADELUS_HAVE_TIME_MONITOR
  using Teuchos::TimeMonitor;
#endif

  using value_type      = typename ZViewType::value_type;
#ifdef PRINT_STATUS
  using execution_space = typename ZViewType::device_type::execution_space;
  using memory_space    = typename ZViewType::device_type::memory_space;
#endif

  double run_secs; // time (in secs) during which the prog ran
  double tsecs;    // intermediate storage of timing info

#ifdef PRINT_STATUS
  printf("Rank %i -- solve_() Begin FwdSolve+BwdSolve+Perm with blksz %d, myrow %d, mycol %d, nprocs_row %d, nprocs_col %d, nrows_matrix %d, ncols_matrix %d, my_rows %d, my_cols %d, my_rhs %d, nrhs %d, value_type %s, execution_space %s, memory_space %s\n", ahandle.get_myrank(), ahandle.get_blksz(), ahandle.get_myrow(), ahandle.get_mycol(), ahandle.get_nprocs_row(), ahandle.get_nprocs_col(), ahandle.get_nrows_matrix(), ahandle.get_ncols_matrix(), ahandle.get_my_rows(), ahandle.get_my_cols(), ahandle.get_my_rhs(), ahandle.get_nrhs(), typeid(value_type).name(), typeid(execution_space).name(), typeid(memory_space).name());
#endif

  {
    tsecs = get_seconds(0.0);

#ifdef ADELUS_PERM_MAT_FORWARD_COPY_TO_HOST
    typename ZViewType::HostMirror h_Z = Kokkos::create_mirror_view( Z );
    typename RHSViewType::HostMirror h_RHS = Kokkos::create_mirror_view( RHS );
    // Bring data to host memory
    Kokkos::deep_copy (h_Z, Z);
    Kokkos::deep_copy (h_RHS, RHS);
#endif

#ifdef ADELUS_HAVE_TIME_MONITOR
    {
      TimeMonitor t(*TimeMonitor::getNewTimer("Adelus: rhs permutation"));
#endif
      // Permute the RHS
#ifdef ADELUS_PERM_MAT_FORWARD_COPY_TO_HOST
      permute_rhs(ahandle, h_RHS, permute);
#else
      permute_rhs(ahandle, RHS, permute);
#endif
#ifdef ADELUS_HAVE_TIME_MONITOR
    }
#endif

#ifdef ADELUS_HAVE_TIME_MONITOR
    {
      TimeMonitor t(*TimeMonitor::getNewTimer("Adelus: forward solve"));
#endif
      //Forward Solve
#ifdef ADELUS_PERM_MAT_FORWARD_COPY_TO_HOST
      forward(ahandle, h_Z, h_RHS);
#else
      forward(ahandle, Z, RHS);
#endif
#ifdef ADELUS_HAVE_TIME_MONITOR
    }
#endif

#ifdef ADELUS_PERM_MAT_FORWARD_COPY_TO_HOST
    // Copy back to device memory
    Kokkos::deep_copy (Z,   h_Z);
    Kokkos::deep_copy (RHS, h_RHS);
#endif

    MPI_Barrier(ahandle.get_comm());

#ifdef ADELUS_HAVE_TIME_MONITOR
    {
      TimeMonitor t(*TimeMonitor::getNewTimer("Adelus: backsolve"));
#endif
      back_solve6(ahandle, Z, RHS);
#ifdef ADELUS_HAVE_TIME_MONITOR
    }
#endif

    MPI_Barrier(ahandle.get_comm());

#ifdef ADELUS_HAVE_TIME_MONITOR
    {
      TimeMonitor t(*TimeMonitor::getNewTimer("Adelus: permutation"));
#endif
      perm1_(ahandle, RHS);
#ifdef ADELUS_HAVE_TIME_MONITOR
    }
#endif

    MPI_Barrier(ahandle.get_comm());

    tsecs = get_seconds(tsecs);

    run_secs = (double) tsecs;

    *secs = run_secs;
    showtime(ahandle.get_comm_id(), ahandle.get_comm(), ahandle.get_myrank(), ahandle.get_nprocs_cube(),
              "Total time in Solve", &run_secs );
  }
}

}//namespace Adelus

#endif
