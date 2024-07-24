// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Details_Profiling.hpp"
#include "Tpetra_Details_Behavior.hpp"
// I don't like including everything, but currently (as of 19 Apr
// 2017), Kokkos does not provide a public header file to get just the
// Profiling hooks.  Users are just supposed to include
// Kokkos_Core.hpp.
#include "Kokkos_Core.hpp"

namespace Tpetra {
namespace Details {

ProfilingRegion::ProfilingRegion ():
  kokkos_region_active_(false)
{

}

ProfilingRegion::ProfilingRegion (const char name[]) {
  kokkos_region_active_ = false;
  if(Behavior::profilingRegionUseKokkosProfiling()){
    kokkos_region_active_ = true;
    ::Kokkos::Profiling::pushRegion(name);
  }
  if(Behavior::profilingRegionUseTeuchosTimers())
    tm = Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(name)));

}

ProfilingRegion::ProfilingRegion (const char name[], const char group[]) {
  kokkos_region_active_ = false;
  const bool timeit = Behavior::timing(group);
  if (timeit)
  {
    if(Behavior::profilingRegionUseKokkosProfiling()){
      kokkos_region_active_ = true;
      ::Kokkos::Profiling::pushRegion(name);
    }
    if(Behavior::profilingRegionUseTeuchosTimers())
      tm = Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(name)));
  }
}

ProfilingRegion::~ProfilingRegion () {
  if(Behavior::profilingRegionUseKokkosProfiling() && kokkos_region_active_){
    ::Kokkos::Profiling::popRegion();
  }
}

} // namespace Details
} // namespace Tpetra


