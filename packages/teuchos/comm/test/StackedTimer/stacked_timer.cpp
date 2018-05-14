// @HEADER
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_StackedTimer.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_DefaultSerialComm.hpp"

#include <thread> // std::this_thread::sleep_for;
Teuchos::RCP<Teuchos::Comm<int> > comm;

int main(int argc, char **argv) {
  int result;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
#ifdef HAVE_MPI
  comm = Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(MPI_COMM_WORLD)));
#else
  comm = Teuchos::rcp(new Teuchos::SerialComm<int>());
#endif
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

}

TEUCHOS_UNIT_TEST(StackedTimer, Basic)
{
  Teuchos::StackedTimer timer("My New Timer");
  timer.start("Total Time");
  {
    for (int i=0; i < 10; ++i) {
    
      timer.start("Assembly");
      std::this_thread::sleep_for(std::chrono::milliseconds{100});   
      timer.stop("Assembly");
      
      timer.start("Solve");
      {
        timer.start("Prec");
        std::this_thread::sleep_for(std::chrono::milliseconds{50});   
        timer.stop("Prec");
        if (rank(*comm) == 0 ) {
          const std::string label = "Rank 0 ONLY";
          timer.start(label);
          std::this_thread::sleep_for(std::chrono::milliseconds{50});
          TEST_ASSERT((timer.findTimer("My New Timer@Solve@Rank 0 ONLY")).running);
          timer.stop(label);
          TEST_ASSERT(not (timer.findTimer("My New Timer@Solve@Rank 0 ONLY")).running);
        } else {
          timer.start("Not Rank 0");
          std::this_thread::sleep_for(std::chrono::milliseconds{50});
          TEST_ASSERT((timer.findTimer("My New Timer@Solve@Not Rank 0")).running);
          timer.stop("Not Rank 0");
          TEST_ASSERT(not (timer.findTimer("My New Timer@Solve@Not Rank 0")).running);
        }
      }
      timer.stop("Solve");
      
    }
  }
  timer.stop("Total Time");
  Teuchos::StackedTimer::OutputOptions options;
  options.output_histogram=true;
  options.num_histogram=3;
  timer.report(std::cout, comm, options);
}


TEUCHOS_UNIT_TEST(StackedTimer, TimeMonitorInteroperability)
{
  const auto precTimer = Teuchos::TimeMonitor::getNewTimer("Prec");
  const auto gmresTimer = Teuchos::TimeMonitor::getNewTimer("GMRES");

  // Test the set and get stacked timer methods on TimeMonitor
  const auto timeMonitorDefaultStackedTimer = Teuchos::TimeMonitor::getStackedTimer();
  const auto timer = Teuchos::rcp(new Teuchos::StackedTimer("StackedTimerTest::TimeMonitorInteroperability"));
  TEST_ASSERT(nonnull(timeMonitorDefaultStackedTimer));
  TEST_ASSERT(nonnull(timer));
  TEST_ASSERT(timeMonitorDefaultStackedTimer != timer);
  Teuchos::TimeMonitor::setStackedTimer(timer);
  TEST_ASSERT(timer == Teuchos::TimeMonitor::getStackedTimer());

  timer->start("Total Time");
  {
    for (int i=0; i < 10; ++i) {
    
      timer->start("Assembly");
      std::this_thread::sleep_for(std::chrono::milliseconds{100});   
      timer->stop("Assembly");
      
      timer->start("Solve");
      {
        {
          Teuchos::TimeMonitor tm(*precTimer);
          std::this_thread::sleep_for(std::chrono::milliseconds{50});
        }
        {
          Teuchos::TimeMonitor tm(*gmresTimer);
          std::this_thread::sleep_for(std::chrono::milliseconds{50});
        }
      }
      timer->stop("Solve");
      
    }
  }
  timer->stop("Total Time");
  
  assert(size(*comm)>0);
  timer->report(std::cout, comm);

#ifdef HAVE_TEUCHOS_ADD_TIME_MONITOR_TO_STACKED_TIMER
  // If time monitor insertion enabled, then test it
#endif
}
