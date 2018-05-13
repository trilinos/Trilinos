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
        if (rank(*comm) != 0 ) {
          timer.start("GMRES");
          std::this_thread::sleep_for(std::chrono::milliseconds{50});
          timer.stop("GMRES");
        } else {
          timer.start("Not GMRES");
          std::this_thread::sleep_for(std::chrono::milliseconds{50});
          timer.stop("Not GMRES");

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
  
  const auto timer = Teuchos::rcp(new Teuchos::StackedTimer("StackedTimerTest::TimeMonitorInteroperability"));

  // Assign timer for TimeMonitor to use
#ifdef HAVE_TEUCHOS_ADD_TIME_MONITOR_TO_STACKED_TIMER
  Teuchos::TimeMonitor::setStackedTimer(timer);
#endif
  
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
}
