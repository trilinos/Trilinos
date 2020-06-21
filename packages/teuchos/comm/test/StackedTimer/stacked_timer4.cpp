// @HEADER
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

#if defined(HAVE_TEUCHOS_KOKKOS_PROFILING) && defined(HAVE_TEUCHOSCORE_KOKKOSCORE)
#include "Kokkos_Core.hpp"
#endif


TEUCHOS_UNIT_TEST(StackedTimer, minmax_hist)
{

  Teuchos::StackedTimer timer("L0");
  timer.stopBaseTimer();

  const Teuchos::RCP<const Teuchos::Comm<int>> comm = Teuchos::DefaultComm<int>::getComm();
  if (comm->getSize() != 4)
    return;
  const int myRank = Teuchos::rank(*comm);

  if ((myRank == 1) || (myRank == 2)) {
    timer.start("T0");
    timer.stop("T0");
    const_cast<Teuchos::BaseTimer*>(timer.findBaseTimer("L0@T0"))->setAccumulatedTime(Teuchos::as<double>(myRank));
  } else {
    timer.start("T1");
    timer.stop("T1");
    const_cast<Teuchos::BaseTimer*>(timer.findBaseTimer("L0@T1"))->setAccumulatedTime(Teuchos::as<double>(myRank));
  }

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
  // Kokkos::finalize_all().
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
#if defined(HAVE_TEUCHOS_KOKKOS_PROFILING) && defined(HAVE_TEUCHOSCORE_KOKKOSCORE)
  Kokkos::initialize(argc,argv);
#endif
  {
    Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
    out.setOutputToRootOnly(0);
  }
  Teuchos::UnitTestRepository::setGloballyReduceTestResult(true);

  auto return_val = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
#if defined(HAVE_TEUCHOS_KOKKOS_PROFILING) && defined(HAVE_TEUCHOSCORE_KOKKOSCORE)
  if (Kokkos::is_initialized())
    Kokkos::finalize_all();
#endif
  return return_val;
}
