// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_PerformanceMonitorBase.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_StackedTimer.hpp"
#include "Teuchos_DefaultComm.hpp"
#include <sstream>
#include <thread> // std::this_thread::sleep_for;
#include <tuple>
#include <regex>
#include <iterator>
#include <limits>

#if defined(HAVE_TEUCHOS_KOKKOS_PROFILING) && defined(HAVE_TEUCHOSCORE_KOKKOS)
#include "Kokkos_Core.hpp"
#endif

// ***************************************************
// These tests must be run on exactly 4 MPI processes!
// ***************************************************

TEUCHOS_UNIT_TEST(StackedTimer, minmax_hist)
{
  Teuchos::StackedTimer timer("L0");
  timer.stopBaseTimer();

  const Teuchos::RCP<const Teuchos::Comm<int>> comm = Teuchos::DefaultComm<int>::getComm();
  TEST_ASSERT(comm->getSize() == 4)
  const int myRank = Teuchos::rank(*comm);

  if ((myRank == 1) || (myRank == 2)) {
    timer.start("T0");
    timer.stop("T0");
    if(myRank == 1) {
      timer.start("T0");
      timer.stop("T0");
      timer.start("T0");
      timer.stop("T0");
    }
    const_cast<Teuchos::BaseTimer*>(timer.findBaseTimer("L0@T0"))->setAccumulatedTime(Teuchos::as<double>(myRank));
  } else {
    timer.start("T1");
    timer.stop("T1");
    if (myRank == 3) {
      timer.start("T1");
      timer.stop("T1");
    }
    const_cast<Teuchos::BaseTimer*>(timer.findBaseTimer("L0@T1"))->setAccumulatedTime(Teuchos::as<double>(myRank));
  }

  // Throws since the mpi aggregation has not been called yet
  TEST_THROW(timer.getMpiAverageTime("L0@T0"),std::runtime_error);
  TEST_THROW(timer.getMpiAverageCount("L0@T1"),std::runtime_error);
  TEST_THROW(timer.isTimer("L0@T1"),std::runtime_error);

  Teuchos::StackedTimer::OutputOptions options;

  out << "\n### Printing default report ###" << std::endl;
  options.output_minmax=true;
  options.output_proc_minmax=true;
  options.output_histogram=true;
  options.num_histogram=2;
  options.output_fraction=false;
  timer.report(out, comm, options);
  {
    std::ostringstream os;
    timer.report(os, comm, options);
    if (myRank == 0) {
      TEST_ASSERT(os.str().find("min=1, max=2, proc min=1, proc max=2") != std::string::npos);
      TEST_ASSERT(os.str().find("<1, 1>") != std::string::npos);
      TEST_ASSERT(os.str().find("min=0, max=3, proc min=0, proc max=3") != std::string::npos);

      constexpr double tol = 10.0*std::numeric_limits<double>::epsilon();
      TEST_FLOATING_EQUALITY(timer.getMpiAverageTime("L0@T0"),1.5,tol);
      TEST_FLOATING_EQUALITY(timer.getMpiAverageTime("L0@T1"),1.5,tol);
      TEST_FLOATING_EQUALITY(timer.getMpiAverageCount("L0@T0"),2.0,tol);
      TEST_FLOATING_EQUALITY(timer.getMpiAverageCount("L0@T1"),1.5,tol);

      TEST_THROW(timer.getMpiAverageTime("INCORRECT TIMER NAME"),std::runtime_error);
      TEST_THROW(timer.getMpiAverageCount("INCORRECT TIMER NAME"),std::runtime_error);

      TEST_ASSERT(timer.isTimer("L0@T0"));
      TEST_ASSERT(timer.isTimer("L0@T1"));
      TEST_ASSERT(!timer.isTimer("INCORRECT TIMER NAME"));
    }
  }
}

// Use our own main to initialize kokkos before calling
// runUnitTestsFromMain(). The kokkos space_time_stack profiler seg
// faults due to inconsistent push/pop of timers in the teuchos unit
// test startup code. By calling initialize here we can use the
// space_time_stack profiler with this unit test.
int main( int argc, char* argv[] )
{
  // Note that the dtor for GlobalMPISession will call
  // Kokkos::finalize().
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
#if defined(HAVE_TEUCHOS_KOKKOS_PROFILING) && defined(HAVE_TEUCHOSCORE_KOKKOS)
  Kokkos::initialize(argc,argv);
#endif
  {
    Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
    out.setOutputToRootOnly(0);
  }
  Teuchos::UnitTestRepository::setGloballyReduceTestResult(true);

  auto return_val = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
#if defined(HAVE_TEUCHOS_KOKKOS_PROFILING) && defined(HAVE_TEUCHOSCORE_KOKKOS)
  if (Kokkos::is_initialized())
    Kokkos::finalize();
#endif
  return return_val;
}
