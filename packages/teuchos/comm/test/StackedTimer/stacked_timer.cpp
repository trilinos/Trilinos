// @HEADER
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_PerformanceMonitorBase.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_StackedTimer.hpp"
#include "Teuchos_DefaultComm.hpp"
#include <sstream>
#include <thread> // std::this_thread::sleep_for;
#include <tuple>
#include <regex>
#include <iterator>

TEUCHOS_UNIT_TEST(PerformanceMonitorBase, UnsortedMergeUnion) {

  const Teuchos::RCP<const Teuchos::Comm<int>> comm = Teuchos::DefaultComm<int>::getComm();

  Teuchos::Array<std::string> a,b, tmp_a, tmp_b;

  a.push_back("foo");
  a.push_back("bar");
  a.push_back("car");

  b.push_back("car");
  b.push_back("bar");
  b.push_back("cat");


  tmp_a=a;
  tmp_b=b;
  Teuchos::unsortedMergePair(tmp_a, tmp_b, Teuchos::Union);
  TEST_EQUALITY(tmp_b.size(),4);
  TEST_EQUALITY(tmp_b[0], "car");
  TEST_EQUALITY(tmp_b[1], "bar");
  TEST_EQUALITY(tmp_b[2], "cat");
  TEST_EQUALITY(tmp_b[3], "foo");
}

TEUCHOS_UNIT_TEST(PerformanceMonitorBase, UnsortedMergeIntersection) {

  const Teuchos::RCP<const Teuchos::Comm<int>> comm = Teuchos::DefaultComm<int>::getComm();

  Teuchos::Array<std::string> a,b, tmp_a, tmp_b;

  a.push_back("foo");
  a.push_back("bar");
  a.push_back("car");

  b.push_back("car");
  b.push_back("bar");
  b.push_back("cat");


  tmp_a=a;
  tmp_b=b;
  Teuchos::unsortedMergePair(tmp_a, tmp_b, Teuchos::Intersection);
  TEST_EQUALITY(tmp_b.size(),2);
  TEST_EQUALITY(tmp_b[0], "car");
  TEST_EQUALITY(tmp_b[1], "bar");
}


TEUCHOS_UNIT_TEST(StackedTimer, Basic)
{
  const Teuchos::RCP<const Teuchos::Comm<int>> comm = Teuchos::DefaultComm<int>::getComm();
  const int myRank = Teuchos::rank(*comm);

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

        // Test different timers on different mpi processes
        if (myRank == 0 ) {
          const std::string label = "Rank 0 ONLY";
          timer.start(label);
          std::this_thread::sleep_for(std::chrono::milliseconds{50});
          TEST_ASSERT((timer.findTimer("Total Time@Solve@Rank 0 ONLY")).running);
          timer.stop(label);
          TEST_ASSERT(not (timer.findTimer("Total Time@Solve@Rank 0 ONLY")).running);
        } else {
          timer.start("Not Rank 0");
          std::this_thread::sleep_for(std::chrono::milliseconds{50});
          TEST_ASSERT((timer.findTimer("Total Time@Solve@Not Rank 0")).running);
          timer.stop("Not Rank 0");
          TEST_ASSERT(not (timer.findTimer("Total Time@Solve@Not Rank 0")).running);
        }
      }
      timer.stop("Solve");

    }
  }
  timer.stop("Total Time");


  TEST_EQUALITY((timer.findTimer("Total Time")).count, 1);
  TEST_EQUALITY((timer.findTimer("Total Time@Assembly")).count, 10);
  TEST_EQUALITY((timer.findTimer("Total Time@Solve")).count, 10);
  TEST_EQUALITY((timer.findTimer("Total Time@Solve@Prec")).count, 10);

  // Test for exception for bad timer name
  TEST_THROW(timer.findTimer("Testing misspelled timer name!"),std::runtime_error);

  // Pre-aggregation
  if (myRank == 0) {
    TEST_EQUALITY((timer.findTimer("Total Time@Solve@Rank 0 ONLY")).count, 10);
  }
  else {
    TEST_EQUALITY((timer.findTimer("Total Time@Solve@Not Rank 0")).count, 10);
  }

  Teuchos::StackedTimer::OutputOptions options;
  options.output_histogram=true;
  options.num_histogram=3;
  options.print_warnings=false;

  // Get the report
  std::stringstream sout1;
  timer.report(sout1, comm, options);

  // Make sure can call report() multiple times, i.e. aggregation
  // resets correctly for each call report()
  std::stringstream sout2;
  timer.report(sout2, comm, options);
  TEST_EQUALITY(sout1.str(),sout2.str());

  // Gold file results (timer name,expected runtime,number of calls)
  std::vector<std::tuple<std::string,double,unsigned long>> lineChecks;
  lineChecks.push_back(std::make_tuple("My New Timer:",0.0,1));
  lineChecks.push_back(std::make_tuple("Total Time:",2.0,1));
  lineChecks.push_back(std::make_tuple("Assembly:",1.0,10));
  lineChecks.push_back(std::make_tuple("Solve:",1.0,10));
  lineChecks.push_back(std::make_tuple("Prec:",0.5,10));

  // Check the report() output. Read the first few lines and parse the
  // expected timer label, the runtime and the counts.
  //
  // * NOTE: The report only combines values to a single MPI process, so
  //         only check on that process.
  // * NOTE: regex not supported in gcc until 4.9. Can drop this check
  //         when Trilinos drops support for gcc 4.8.
#if !defined(__GNUC__) \
    || ( defined(__GNUC__) && (__GNUC__ > 4) ) \
    || ( defined(__GNUC__) && (__GNUC__ == 4) && (__GNUC__MINOR__ > 8) )

  if (myRank == 0) {
    const double timerTolerance = 0.25; // +- 0.25 seconds
    std::istringstream is(sout1.str());
    for (const auto& check : lineChecks) {

      std::string line;
      std::getline(is,line);
      std::smatch regexSMatch;
      std::regex timerName(std::get<0>(check));
      std::regex_search(line,regexSMatch,timerName);
      TEST_ASSERT(not regexSMatch.empty());

      // Split string to get time and count
      std::regex delimiter(":\\s|\\s\\[|\\]\\s");
      std::sregex_token_iterator tok(line.begin(), line.end(),delimiter,-1);

      const std::string timeAsString = (++tok)->str();
      const double time = std::stod(timeAsString);
      TEST_FLOATING_EQUALITY(time,std::get<1>(check),timerTolerance);

      const std::string countAsString = (++tok)->str();
      const unsigned long count = std::stoul(countAsString);
      TEST_EQUALITY(count,std::get<2>(check));
    }
  }
#endif

  // Print to screen
  out << "\n### Printing default report ###" << std::endl;
  Teuchos::StackedTimer::OutputOptions defaultOptions;
  timer.report(out, comm, defaultOptions);

  // Test some options
  out << "\n### Printing aligned_column with timers names on left ###" << std::endl;
  options.output_fraction = true;
  options.output_total_updates = true;
  options.output_minmax = true;
  options.output_histogram = true;
  options.num_histogram = 3;
  options.align_columns = true;
  timer.report(out, comm, options);

  // Toggle names before values
  TEST_EQUALITY(options.print_names_before_values,true);
  out << "\n### Printing aligned_column with timers names on right ###" << std::endl;
  options.print_names_before_values = false;
  timer.report(out, comm, options);
}


TEUCHOS_UNIT_TEST(StackedTimer, TimeMonitorInteroperability)
{
  const Teuchos::RCP<const Teuchos::Comm<int>> comm = Teuchos::DefaultComm<int>::getComm();

  const auto diffTimer = Teuchos::TimeMonitor::getNewTimer("Diffusion Term");
  const auto rxnTimer = Teuchos::TimeMonitor::getNewTimer("Reaction Term");
  const auto precTimer = Teuchos::TimeMonitor::getNewTimer("Prec");
  const auto gmresTimer = Teuchos::TimeMonitor::getNewTimer("GMRES");

  // Test the set and get stacked timer methods on TimeMonitor
  const auto timeMonitorDefaultStackedTimer = Teuchos::TimeMonitor::getStackedTimer();
  const auto timer = Teuchos::rcp(new Teuchos::StackedTimer("TM:Interoperability"));
  TEST_ASSERT(nonnull(timeMonitorDefaultStackedTimer));
  TEST_ASSERT(nonnull(timer));
  TEST_ASSERT(timeMonitorDefaultStackedTimer != timer);
  Teuchos::TimeMonitor::setStackedTimer(timer);
  TEST_ASSERT(timer == Teuchos::TimeMonitor::getStackedTimer());

  timer->start("Total Time");
  {
    for (int i=0; i < 10; ++i) {

      timer->start("Assembly");
      {
        {
          Teuchos::TimeMonitor tm(*diffTimer);
          std::this_thread::sleep_for(std::chrono::milliseconds{25});
        }
        {
          Teuchos::TimeMonitor tm(*rxnTimer);
          std::this_thread::sleep_for(std::chrono::milliseconds{75});
        }
        // Remainder
        std::this_thread::sleep_for(std::chrono::milliseconds{100});
      }
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
        // Remainder
        std::this_thread::sleep_for(std::chrono::milliseconds{100});
      }
      timer->stop("Solve");
      std::this_thread::sleep_for(std::chrono::milliseconds{100});
    }
  }
  timer->stop("Total Time");

  assert(size(*comm)>0);

  TEST_EQUALITY((timer->findTimer("Total Time")).count, 1);
  TEST_EQUALITY((timer->findTimer("Total Time@Assembly")).count, 10);

  // Make sure the TimeMonitor added the timers
#ifdef HAVE_TEUCHOS_ADD_TIME_MONITOR_TO_STACKED_TIMER
  TEST_EQUALITY((timer->findTimer("Total Time@Solve@Prec")).count, 10);
  TEST_EQUALITY((timer->findTimer("Total Time@Solve@GMRES")).count, 10);
#endif

  Teuchos::StackedTimer::OutputOptions options;
  out << "\n### Printing default report ###" << std::endl;
  options.output_histogram=true;
  options.num_histogram=3;
  options.output_fraction=true;
  timer->report(out, comm, options);

  out << "\n### Printing aligned_column with timers names on left ###" << std::endl;
  options.align_columns = true;
  timer->report(out, comm, options);

  out << "\n### Printing aligned_column with timers names on right ###" << std::endl;
  // options.print_names_before_values=false requires that
  // options.align_output=true. The code will automatically fix this
  // and print a warning if warnings are enabled. Testing this here by
  // specifying the incorrect logic.
  options.align_columns = false;
  options.print_names_before_values = false;
  timer->report(out, comm, options);

  //Testing limited number of levels in printing
  out << "\n### Printing with max_levels=2 ###" << std::endl;
  options.max_levels=2;
  options.align_columns = true;
  options.print_names_before_values = true;
  timer->report(out, comm, options);
}

// Overlapping timers are not allowed in a StackedTimer, but are in
// TimeMonitor. Since StackedTimer is automatically used in
// TimeMonitor by default, we have seen this error - a throw from the
// stacked timer. In every instance so far, the intention was not to
// actually overlap but a constructor/destructor ordering issue
// (usually involving RCPs). To prevent tests from failing,
// StackedTimer now automatically shuts itself off if it detects
// overlaped timers in a TimeMonitor instance, reports a warning on
// how to fix and allows the code to continue runnning. Where this has
// occurred in Trilinos is when a TimeMonitor object is stored in an
// RCP and then the RCP is reassigned to a new timer. The intention
// was to stop one and start another. But the destruction of one and
// the creation of the new one occurs in the wrong order. This test
// demonstrates the issue.
TEUCHOS_UNIT_TEST(StackedTimer, OverlappingTimersException)
{
  Teuchos::StackedTimer timer("My Timer");
  timer.start("Outer");
  timer.start("Inner");
  // Should stop inner before outer
  TEST_THROW(timer.stop("Outer"),std::runtime_error);
}


#ifdef HAVE_TEUCHOS_ADD_TIME_MONITOR_TO_STACKED_TIMER
TEUCHOS_UNIT_TEST(StackedTimer, OverlappingTimersViaRCP)
{
  const auto precTimer = Teuchos::TimeMonitor::getNewTimer("Prec");
  const auto gmresTimer = Teuchos::TimeMonitor::getNewTimer("GMRES");

  Teuchos::RCP<Teuchos::TimeMonitor> timer = Teuchos::rcp(new Teuchos::TimeMonitor(*precTimer));
  timer = Teuchos::rcp(new Teuchos::TimeMonitor(*gmresTimer));

  TEST_ASSERT(is_null(Teuchos::TimeMonitor::getStackedTimer()));
}
#endif


