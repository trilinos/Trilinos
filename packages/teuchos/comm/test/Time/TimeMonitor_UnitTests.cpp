// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_GlobalMPISession.hpp"

// slowLoop does not reliably make a timer nonzero (RHEL6, gcc 4.4.7, OpenMPI 1.5.4).
// Thus, I'm introducing headers to make sleep() available.
#ifndef ICL
#include <unistd.h>
#else
void sleep(int sec)
{
  Sleep(sec);
}
#endif

#ifdef _MSC_VER
#pragma comment(lib, "Ws2_32.lib")
# include <Winsock2.h>
# include <process.h>
void sleep(int sec)
{
  Sleep(sec * 1000);
}
#endif

#ifdef HAVE_TEUCHOS_TIMER_KOKKOS_FENCE
#include "Kokkos_Core.hpp"
#endif

namespace {

  void func_time_monitor1()
  {
    TEUCHOS_FUNC_TIME_MONITOR("FUNC_TIME_MONITOR1");
  }


  void func_time_monitor2()
  {
    TEUCHOS_FUNC_TIME_MONITOR("FUNC_TIME_MONITOR2");
    {
      TEUCHOS_FUNC_TIME_MONITOR_DIFF("FUNC_TIME_MONITOR2_inner", inner);
    }
  }

  // This function is commented out because no code in this file
  // actually uses it.  Commenting it out avoids compiler warnings
  // ("unused function").
#if 0
  // Slowly compute the n-th Fibonacci number.  This gives timers
  // something to time.  Be careful not to make n too large, or you'll
  // run out of stack space.
  int
  fib (const int n)
  {
    if (n <= 0) {
      return 0;
    } else if (n == 1) {
      return 1;
    }
    else {
      // You should never compute the n-th Fibonacci number like this.
      // This is exponentially slow in n.  The point of using a slow
      // algorithm like this is to exercise timers.
      return fib (n-1) + fib (n-2);
    }
  }
#endif // 0

  // Do a number of arithmetic operations proportional to n^3, in
  // order to have something to time.  Unlike the recursive function
  // fib() commented out above, this function shouldn't take up a lot
  // of stack space.
  double
  slowLoop (const size_t n)
  {
    const double inc = 1 / static_cast<double> (n);
    double x = 1.0;

    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j) {
        for (size_t k = 0; k < n; ++k) {
          x = x + inc;
        }
      }
    }
    return x;
  }

} // namespace (anonymous)


namespace Teuchos {
  //
  // Basic TimeMonitor test: create and exercise a timer on all (MPI)
  // processes, and make sure that TimeMonitor::summarize() reports it.
  //
  TEUCHOS_UNIT_TEST( TimeMonitor, FUNC_TIME_MONITOR  )
  {
    using Teuchos::ParameterList;
    using Teuchos::parameterList;
    using Teuchos::RCP;

    func_time_monitor1 (); // Function to time.

    { // Repeat test for default output format.
      std::ostringstream oss;
      TimeMonitor::summarize (oss);

      // Echo summarize() output to the FancyOStream out (which is a
      // standard unit test argument).  Output should only appear in
      // show-all-test-details mode.
      out << oss.str () << std::endl;

      // Make sure that the timer's name shows up in the output.
      if (Teuchos::GlobalMPISession::getRank() == 0) {
        const size_t substr_i = oss.str ().find ("FUNC_TIME_MONITOR1");
        TEST_INEQUALITY(substr_i, std::string::npos);
      }
    }

    { // Repeat test for YAML output, compact style.
      std::ostringstream oss;
      RCP<ParameterList> reportParams =
        parameterList (* (TimeMonitor::getValidReportParameters ()));
      reportParams->set ("Report format", "YAML");
      reportParams->set ("YAML style", "compact");
      TimeMonitor::report (oss, reportParams);

      // Echo output to the FancyOStream out (which is a standard unit
      // test argument).  Output should only appear in "show all test
      // details" mode.
      out << oss.str () << std::endl;

      // Make sure that the timer's name shows up in the output.
      const size_t substr_i = oss.str ().find ("FUNC_TIME_MONITOR1");
      TEST_INEQUALITY(substr_i, std::string::npos);
    }

    { // Repeat test for YAML output, spacious style.
      std::ostringstream oss;
      RCP<ParameterList> reportParams =
        parameterList (* (TimeMonitor::getValidReportParameters ()));
      reportParams->set ("Report format", "YAML");
      reportParams->set ("YAML style", "spacious");
      TimeMonitor::report (oss, reportParams);

      // Echo output to the FancyOStream out (which is a standard unit
      // test argument).  Output should only appear in "show all test
      // details" mode.
      out << oss.str () << std::endl;

      // Make sure that the timer's name shows up in the output.
      const size_t substr_i = oss.str ().find ("FUNC_TIME_MONITOR1");
      TEST_INEQUALITY(substr_i, std::string::npos);
    }

    // This sets up for the next unit test.
    TimeMonitor::clearCounters ();
  }

  //
  // Test for TimeMonitor's enableTimer and disableTimer methods.
  //
  TEUCHOS_UNIT_TEST( TimeMonitor, enableTimer )
  {
    using Teuchos::Array;
    using Teuchos::OSTab;
    using Teuchos::ParameterList;
    using Teuchos::parameterList;
    using Teuchos::RCP;
    using Teuchos::rcpFromRef;
    using Teuchos::Time;
    using Teuchos::TimeMonitor;
    using std::endl;
    typedef Teuchos::Array<RCP<Time> >::size_type size_type;

    out << "Testing TimeMonitor's disableTimer() and enableTimer() methods"
        << endl;
    OSTab (rcpFromRef (out));

    out << "Creating timers" << endl;
    const int numTrials = 5;
    const int numTimers = 3;
    Array<RCP<Time> > timers (numTimers);
    for (size_type i = 0; i < numTimers; ++i) {
      std::ostringstream os; // construct timer name
      os << "Timer " << i;
      timers[i] = TimeMonitor::getNewTimer (os.str ());
    }

    out << "Running all timers without disabling any" << endl;
    // The actual number of operations in slowloop is proportional to
    // the cube of the loop length.  Adjust loopLength as necessary to
    // ensure the timer reports a nonzero elapsed time for each of the
    // invocations.
    const size_t loopLength = 25;
    for (int k = 0; k < numTrials; ++k) {
      for (size_type i = 0; i < numTimers; ++i) {
        TimeMonitor timeMon (* timers[i]);
        slowLoop (loopLength);
      }
    }
    for (size_type i = 0; i < numTimers; ++i) {
      TEST_EQUALITY( timers[i]->numCalls(), numTrials );
    }

    out << "Disabling one timer and trying again" << endl;
    // Disable timers[0] only, and repeat the above loops.
    TEST_NOTHROW( TimeMonitor::disableTimer ("Timer 0") );
    for (int k = 0; k < numTrials; ++k) {
      for (size_type i = 0; i < numTimers; ++i) {
        TimeMonitor timeMon (* timers[i]);
        slowLoop (loopLength);
      }
    }
    TEST_EQUALITY( timers[0]->numCalls(), numTrials );
    for (size_type i = 1; i < numTimers; ++i) {
      TEST_EQUALITY( timers[i]->numCalls(), 2*numTrials );
    }

    out << "Reenabling the timer and trying again" << endl;
    // Enable timers[0] and repeat the above loops.
    TEST_NOTHROW( TimeMonitor::enableTimer ("Timer 0") );
    for (int k = 0; k < numTrials; ++k) {
      for (size_type i = 0; i < numTimers; ++i) {
        TimeMonitor timeMon (* timers[i]);
        slowLoop (loopLength);
      }
    }
    TEST_EQUALITY( timers[0]->numCalls(), 2*numTrials );
    for (size_type i = 1; i < numTimers; ++i) {
      TEST_EQUALITY( timers[i]->numCalls(), 3*numTrials );
    }

    out << "Test that summarize() reports enabled and disabled timers" << endl;
    // Make sure that summarize() reports all timers.  Disabling a
    // timer must _not_ exclude it from the list of timers printed by
    // summarize().  Disable a different timer this time just for fun.
    TEST_NOTHROW( TimeMonitor::disableTimer ("Timer 1") );
    {
      std::ostringstream oss;
      TimeMonitor::summarize (oss);

      // Echo summarize() output to the FancyOStream out (which is a
      // standard unit test argument).  Output should only appear in
      // show-all-test-details mode.
      out << oss.str () << std::endl;

      // Make sure that each timer's name shows up in the output.
      for (size_type i = 0; i < numTimers; ++i) {
        const std::string name = timers[i]->name ();
        if (Teuchos::GlobalMPISession::getRank() == 0) {
          const size_t substr_i = oss.str ().find (name);
          TEST_INEQUALITY(substr_i, std::string::npos);
        }
      }
    }

    // This sets up for the next unit test.
    TimeMonitor::clearCounters ();
  }


  //
  // Test correct quoting of labels for TimeMonitor's YAML output.
  //
  TEUCHOS_UNIT_TEST( TimeMonitor, YamlLabelQuoting )
  {
    using Teuchos::Array;
    using Teuchos::ParameterList;
    using Teuchos::parameterList;
    using Teuchos::RCP;
    using Teuchos::Time;
    typedef Array<std::string>::size_type size_type;

    Array<std::string> inputLabels, outputLabels;

    // Make sure to exercise things that don't need quoting, like
    // spaces and certain punctuation, as well as things that do need
    // quoting, like colons, inner double quotes, and backslashes.
    inputLabels.push_back ("NoQuotingNeeded");
    inputLabels.push_back ("No quoting needed");
    inputLabels.push_back ("\"AlreadyQuotedNoQuotingNeeded\"");
    inputLabels.push_back ("\"Already quoted, no quoting needed\"");
    inputLabels.push_back ("\"Already quoted: quoting needed\"");
    inputLabels.push_back ("NotQuoted:QuotingNeeded");
    inputLabels.push_back ("Not quoted: quoting needed");
    // Test both individual double quotes, and pairs of double quotes.
    inputLabels.push_back ("Not quoted \" quoting needed");
    inputLabels.push_back ("Not quoted \" \" quoting needed");
    inputLabels.push_back ("\"Already quoted \" quoting needed\"");
    inputLabels.push_back ("\"Already quoted \" \" quoting needed\"");
    // Remember that in C strings, a double backslash turns into a
    // single backslash.  Our YAML output routine should turn each
    // single backslash back into a double backslash.
    inputLabels.push_back ("Not quoted \\ quoting needed");
    inputLabels.push_back ("Not quoted \\\\ quoting needed");
    inputLabels.push_back ("Not quoted \\ \\ quoting needed");
    inputLabels.push_back ("\"Already quoted \\ quoting needed\"");
    inputLabels.push_back ("\"Already quoted \\\\ quoting needed\"");
    inputLabels.push_back ("\"Already quoted \\ \\ quoting needed\"");

    outputLabels.push_back ("NoQuotingNeeded");
    outputLabels.push_back ("No quoting needed");
    outputLabels.push_back ("\"AlreadyQuotedNoQuotingNeeded\"");
    outputLabels.push_back ("\"Already quoted, no quoting needed\"");
    outputLabels.push_back ("\"Already quoted: quoting needed\"");
    outputLabels.push_back ("\"NotQuoted:QuotingNeeded\"");
    outputLabels.push_back ("\"Not quoted: quoting needed\"");
    outputLabels.push_back ("\"Not quoted \\\" quoting needed\"");
    outputLabels.push_back ("\"Not quoted \\\" \\\" quoting needed\"");
    outputLabels.push_back ("\"Already quoted \\\" quoting needed\"");
    outputLabels.push_back ("\"Already quoted \\\" \\\" quoting needed\"");
    outputLabels.push_back ("\"Not quoted \\\\ quoting needed\"");
    outputLabels.push_back ("\"Not quoted \\\\\\\\ quoting needed\"");
    outputLabels.push_back ("\"Not quoted \\\\ \\\\ quoting needed\"");
    outputLabels.push_back ("\"Already quoted \\\\ quoting needed\"");
    outputLabels.push_back ("\"Already quoted \\\\\\\\ quoting needed\"");
    outputLabels.push_back ("\"Already quoted \\\\ \\\\ quoting needed\"");

    // Sanity check.
    TEUCHOS_TEST_FOR_EXCEPTION(
      inputLabels.size () != outputLabels.size (),
      std::logic_error,
      "The number of input labels is different than the number of output labels."
      "  Please ask a Teuchos developer to make sure that every test input "
      "label has a corresponding output label.");

    Array<RCP<Time> > timers;
    for (size_type i = 0; i < inputLabels.size (); ++i) {
      timers.push_back (TimeMonitor::getNewCounter (inputLabels[i]));
    }

    // The actual number of operations in the loop is proportional to
    // the cube of the loop length.  Adjust the quantities below as
    // necessary to ensure the timer reports a nonzero elapsed time
    // for each of the invocations.
    const size_t loopLength = 25;
    for (int k = 0; k < 3; ++k) {
      for (size_type i = 0; i < timers.size (); ++i) {
        TimeMonitor timeMon (* timers[i]);
        slowLoop (loopLength);
      }
    }

    { // YAML output, compact style.
      std::ostringstream oss;
      RCP<ParameterList> reportParams =
        parameterList (* (TimeMonitor::getValidReportParameters ()));
      reportParams->set ("Report format", "YAML");
      reportParams->set ("YAML style", "compact");
      TimeMonitor::report (oss, reportParams);

      // Echo output to the FancyOStream out (which is a standard unit
      // test argument).  Output should only appear in "show all test
      // details" mode.
      out << oss.str () << std::endl;

      // Make sure that all timer labels appear correctly in the output.
      for (size_type i = 0; i < inputLabels.size(); ++i) {
        const size_t pos = oss.str ().find (outputLabels[i]);
        TEST_INEQUALITY(pos, std::string::npos);
      }
    }

    { // YAML output, spacious style.
      std::ostringstream oss;
      RCP<ParameterList> reportParams =
        parameterList (* (TimeMonitor::getValidReportParameters ()));
      reportParams->set ("Report format", "YAML");
      reportParams->set ("YAML style", "spacious");
      TimeMonitor::report (oss, reportParams);

      // Echo output to the FancyOStream out (which is a standard unit
      // test argument).  Output should only appear in "show all test
      // details" mode.
      out << oss.str () << std::endl;

      // Make sure that all timer labels appear correctly in the output.
      for (size_type i = 0; i < inputLabels.size(); ++i) {
        const size_t pos = oss.str ().find (outputLabels[i]);
        TEST_INEQUALITY(pos, std::string::npos);
      }
    }

    // This sets up for the next unit test.
    TimeMonitor::clearCounters ();
  }


  //
  // Test filtering of timer labels.
  //
  TEUCHOS_UNIT_TEST( TimeMonitor, TimerLabelFiltering )
  {
    using Teuchos::Array;
    using Teuchos::ParameterList;
    using Teuchos::parameterList;
    using Teuchos::RCP;
    using Teuchos::Time;
    typedef Array<std::string>::size_type size_type;

    // Filters to use in the test.
    Array<std::string> filters;
    filters.push_back ("Foo:");
    filters.push_back ("Bar:");
    filters.push_back ("Baz:");

    // All the timer labels.
    Array<std::string> labels;
    labels.push_back ("Foo: timer 1");
    labels.push_back ("Foo: timer 2");
    labels.push_back ("Foo: timer 3");
    labels.push_back ("Bar: timer 1");
    labels.push_back ("Bar: timer 2");
    labels.push_back ("Baz: timer 1");
    labels.push_back ("Xyzzy");
    labels.push_back ("This is not a pipe");
    labels.push_back ("You should not see this");

    Array<Array<std::string> > outLabels (3);
    // Label(s) that should be printed for filters[0]
    outLabels[0].push_back ("Foo: timer 1");
    outLabels[0].push_back ("Foo: timer 2");
    outLabels[0].push_back ("Foo: timer 3");
    // Label(s) that should be printed for filters[1]
    outLabels[1].push_back ("Bar: timer 1");
    outLabels[1].push_back ("Bar: timer 2");
    // Label(s) that should be printed for filters[2]
    outLabels[2].push_back ("Baz: timer 1");

    // Labels that should not be printed for any of the filters below.
    Array<std::string> otherLabels;
    otherLabels.push_back ("Xyzzy");
    otherLabels.push_back ("This is not a pipe");
    otherLabels.push_back ("You should not see this");

    Array<RCP<Time> > timers;
    for (size_type i = 0; i < labels.size (); ++i) {
      timers.push_back (TimeMonitor::getNewCounter (labels[i]));
    }

    // The actual number of operations in the loop is proportional to
    // the cube of the loop length.  Adjust the quantities below as
    // necessary to ensure the timer reports a nonzero elapsed time
    // for each of the invocations.
    const size_t loopLength = 25;
    for (int k = 0; k < 3; ++k) {
      for (size_type i = 0; i < timers.size (); ++i) {
        TimeMonitor timeMon (* timers[i]);
        slowLoop (loopLength);
      }
    }

    try {
      // FIXME (mfh 21 Aug 2012) We don't yet have a test ensuring that
      // the filter only selects at the beginning of the timer label.

      // Test for each filter.
      for (size_type i = 0; i < filters.size (); ++i) {
        { // Default (tabular) output format.
          std::ostringstream oss;
          RCP<ParameterList> reportParams =
            parameterList (* (TimeMonitor::getValidReportParameters ()));
          TimeMonitor::report (oss, filters[i], reportParams);

          // Echo output to the FancyOStream out (which is a standard unit
          // test argument).  Output should only appear in "show all test
          // details" mode.
          out << oss.str () << std::endl;

          // Check whether the labels that were supposed to be printed
          // were actually printed.
          if (Teuchos::GlobalMPISession::getRank() == 0) {
            for (size_type j = 0; j < outLabels[i].size(); ++j) {
              const size_t pos = oss.str ().find (outLabels[i][j]);
              TEST_INEQUALITY(pos, std::string::npos);
            }
          }

          // Check whether the labels that were _not_ supposed to be
          // printed were actually printed.
          //
          // First, check the labels that should only be printed with
          // the other filters.
          for (size_type ii = 0; ii < outLabels.size(); ++ii) {
            if (ii != i) {
              for (size_type j = 0; j < outLabels[ii].size(); ++j) {
                const size_t pos = oss.str ().find (outLabels[ii][j]);
                TEST_EQUALITY(pos, std::string::npos);
              }
            }
          }
          // Next, check the labels that should not be printed for any
          // filters.
          for (size_type j = 0; j < otherLabels.size(); ++j) {
            const size_t pos = oss.str ().find (otherLabels[j]);
            TEST_EQUALITY(pos, std::string::npos);
          }
        }

        { // YAML output, compact style.
          std::ostringstream oss;
          RCP<ParameterList> reportParams =
            parameterList (* (TimeMonitor::getValidReportParameters ()));
          reportParams->set ("Report format", "YAML");
          reportParams->set ("YAML style", "compact");
          TimeMonitor::report (oss, filters[i], reportParams);

          // Echo output to the FancyOStream out (which is a standard unit
          // test argument).  Output should only appear in "show all test
          // details" mode.
          out << oss.str () << std::endl;

          // Check whether the labels that were supposed to be printed
          // were actually printed.
          for (size_type j = 0; j < outLabels[i].size(); ++j) {
            const size_t pos = oss.str ().find (outLabels[i][j]);
            TEST_INEQUALITY(pos, std::string::npos);
          }

          // Check whether the labels that were _not_ supposed to be
          // printed were actually printed.
          //
          // First, check the labels that should only be printed with
          // the other filters.
          for (size_type ii = 0; ii < outLabels.size(); ++ii) {
            if (ii != i) {
              for (size_type j = 0; j < outLabels[ii].size(); ++j) {
                const size_t pos = oss.str ().find (outLabels[ii][j]);
                TEST_EQUALITY(pos, std::string::npos);
              }
            }
          }
          // Next, check the labels that should not be printed for any
          // filters.
          for (size_type j = 0; j < otherLabels.size(); ++j) {
            const size_t pos = oss.str ().find (otherLabels[j]);
            TEST_EQUALITY(pos, std::string::npos);
          }
        }

        { // YAML output, spacious style.
          std::ostringstream oss;
          RCP<ParameterList> reportParams =
            parameterList (* (TimeMonitor::getValidReportParameters ()));
          reportParams->set ("Report format", "YAML");
          reportParams->set ("YAML style", "spacious");
          TimeMonitor::report (oss, filters[i], reportParams);

          // Echo output to the FancyOStream out (which is a standard unit
          // test argument).  Output should only appear in "show all test
          // details" mode.
          out << oss.str () << std::endl;

          // Check whether the labels that were supposed to be printed
          // were actually printed.
          for (size_type j = 0; j < outLabels[i].size(); ++j) {
            const size_t pos = oss.str ().find (outLabels[i][j]);
            TEST_INEQUALITY(pos, std::string::npos);
          }

          // Check whether the labels that were _not_ supposed to be
          // printed were actually printed.
          //
          // First, check the labels that should only be printed with
          // the other filters.
          for (size_type ii = 0; ii < outLabels.size(); ++ii) {
            if (ii != i) {
              for (size_type j = 0; j < outLabels[ii].size(); ++j) {
                const size_t pos = oss.str ().find (outLabels[ii][j]);
                TEST_EQUALITY(pos, std::string::npos);
              }
            }
          }
          // Next, check the labels that should not be printed for any
          // filters.
          for (size_type j = 0; j < otherLabels.size(); ++j) {
            const size_t pos = oss.str ().find (otherLabels[j]);
            TEST_EQUALITY(pos, std::string::npos);
          }
        }
      }
    }
    catch (...) {
      // Make sure to clear the counters, so that they don't pollute
      // the remaining tests.  (The Teuchos unit test framework may
      // catch any exceptions that the above code throws, but allow
      // the remaining tests to continue.)
      TimeMonitor::clearCounters ();
      throw;
    }

    // This sets up for the next unit test.
    TimeMonitor::clearCounters ();
  }



  //
  // TimeMonitor nested timers test: create two timers on all (MPI)
  // processes, use the second inside the scope of the first, and make
  // sure that TimeMonitor::summarize() reports both timers.
  //
  TEUCHOS_UNIT_TEST( TimeMonitor, FUNC_TIME_MONITOR_tested  )
  {
    using Teuchos::ParameterList;
    using Teuchos::parameterList;
    using Teuchos::RCP;

    func_time_monitor2 ();

    {
      std::ostringstream oss;
      TimeMonitor::summarize (oss);

      // Echo summarize() output to the FancyOStream out (which is a
      // standard unit test argument).  Output should only appear in
      // show-all-test-details mode.
      out << oss.str() << std::endl;

      if (Teuchos::GlobalMPISession::getRank() == 0) {
        const size_t substr_i = oss.str().find ("FUNC_TIME_MONITOR2");
        TEST_INEQUALITY(substr_i, std::string::npos);
        const size_t substr_inner_i = oss.str().find ("FUNC_TIME_MONITOR2_inner");
        TEST_INEQUALITY(substr_inner_i, std::string::npos);
      }
    }

    { // Repeat test for YAML output, compact style.
      std::ostringstream oss;
      RCP<ParameterList> reportParams =
        parameterList (* (TimeMonitor::getValidReportParameters ()));
      reportParams->set ("Report format", "YAML");
      reportParams->set ("YAML style", "compact");
      TimeMonitor::report (oss, reportParams);

      // Echo output to the FancyOStream out (which is a standard unit
      // test argument).  Output should only appear in "show all test
      // details" mode.
      out << oss.str () << std::endl;

      const size_t substr_i = oss.str().find ("FUNC_TIME_MONITOR2");
      TEST_INEQUALITY(substr_i, std::string::npos);
      const size_t substr_inner_i = oss.str().find ("FUNC_TIME_MONITOR2_inner");
      TEST_INEQUALITY(substr_inner_i, std::string::npos);
    }

    { // Repeat test for YAML output, spacious style.
      std::ostringstream oss;
      RCP<ParameterList> reportParams =
        parameterList (* (TimeMonitor::getValidReportParameters ()));
      reportParams->set ("Report format", "YAML");
      reportParams->set ("YAML style", "spacious");
      TimeMonitor::report (oss, reportParams);

      // Echo output to the FancyOStream out (which is a standard unit
      // test argument).  Output should only appear in "show all test
      // details" mode.
      out << oss.str () << std::endl;

      const size_t substr_i = oss.str().find ("FUNC_TIME_MONITOR2");
      TEST_INEQUALITY(substr_i, std::string::npos);
      const size_t substr_inner_i = oss.str().find ("FUNC_TIME_MONITOR2_inner");
      TEST_INEQUALITY(substr_inner_i, std::string::npos);
    }

    // This sets up for the next unit test.
    TimeMonitor::clearCounters ();
  }

  //
  // Test whether TimeMonitor::summarize() does the right thing when
  // different MPI processes create different sets of timers.  This
  // test is only meaningful if running with more than one MPI
  // process, but it should also pass if running in a non-MPI build or
  // with only one MPI process.
  //
  TEUCHOS_UNIT_TEST( TimeMonitor, SUMMARIZE_diffTimerSets )
  {
    RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm ();
    const int numProcs = comm->getSize ();
    const int myRank = comm->getRank ();

    // If MPI has not been initialized, turn the MPI communicator into
    // a "serial" communicator.  This may not be necessary when using
    // the Teuchos Unit Test Framework.
#ifdef HAVE_MPI
    {
      int mpiHasBeenInitialized = 0;
      MPI_Initialized (&mpiHasBeenInitialized);
      if (! mpiHasBeenInitialized) {
        comm = Teuchos::DefaultComm<int>::getDefaultSerialComm (null);
      }
    }
#endif // HAVE_MPI

    // Store the output of TimeMonitor::summarize() here.
    std::ostringstream oss;

    // Timer A gets created on all MPI processes.
    // Timer B only gets created on Proc 0.
    RCP<Time> timerA = TimeMonitor::getNewCounter ("Timer A");
    RCP<Time> timerB;
    if (myRank == 0) {
      timerB = TimeMonitor::getNewCounter ("Timer B");
    }
    else {
      timerB = null; // True anyway, but I want to make this explicit.
    }

    // The actual number of operations in the loop is proportional to
    // the cube of the loop length.  Adjust the quantities below as
    // necessary to ensure the timer reports a nonzero elapsed time
    // for each of the invocations.
    const size_t timerA_loopLength = 150;
    const size_t timerB_loopLength = 200;

    // Timer A gets a call count of 3.
    for (size_t k = 0; k < 3; ++k) {
      TimeMonitor monitorA (*timerA);
      slowLoop (timerA_loopLength);
    }
    if (myRank == 0) { // Timer B gets a call count of 1 on Proc 0.
      TimeMonitor monitorB (*timerB);
      slowLoop (timerB_loopLength);
    }

    const bool alwaysWriteLocal = false; // the default
    const bool writeGlobalStats = true;  // the default
    const bool writeZeroTimers = true;  // the default
    TimeMonitor::summarize (oss, alwaysWriteLocal, writeGlobalStats,
                            writeZeroTimers, Intersection);

    // Echo summarize() output to the FancyOStream out (which is a
    // standard unit test argument).  Output should only appear in
    // show-all-test-details mode.
    out << std::endl << "Testing intersection of timers:" << std::endl
        << oss.str() << std::endl;

    // Since setOp == Intersection, only Timer A should be reported,
    // unless there is only one (MPI) process.
    if (Teuchos::GlobalMPISession::getRank() == 0) {
     size_t substr_i = oss.str().find ("Timer A");
     TEST_INEQUALITY(substr_i, std::string::npos);
     if (numProcs > 1) {
       substr_i = oss.str().find ("Timer B");
       TEST_EQUALITY(substr_i, std::string::npos);
     }
    }

    // Now call summarize(), but ask for a union of timers.
    std::ostringstream ossUnion;
    TimeMonitor::summarize (ossUnion, alwaysWriteLocal, writeGlobalStats,
                            writeZeroTimers, Union);

    // Echo summarize() output to the FancyOStream out (which is a
    // standard unit test argument).  Output should only appear in
    // show-all-test-details mode.
    out << std::endl << "Testing union of timers:" << std::endl
        << ossUnion.str() << std::endl;

    // Since setOp == Union, both Timer A and Timer B should be
    // reported.
    if (Teuchos::GlobalMPISession::getRank() == 0) {
      size_t substr_i = ossUnion.str().find ("Timer A");
      TEST_INEQUALITY(substr_i, std::string::npos);
      substr_i = ossUnion.str().find ("Timer B");
      TEST_INEQUALITY(substr_i, std::string::npos);
    }

    // This sets up for the next unit test.
    TimeMonitor::clearCounters ();
  }

  //
  // Test whether the option to filter out zero timers works, for the
  // case that all (MPI) processes have the same call counts.
  //
  TEUCHOS_UNIT_TEST( TimeMonitor, FILTER_ZERO_TIMERS_sameParallelCallCounts )
  {
    // Store the output of TimeMonitor::summarize() here.
    std::ostringstream oss;

    RCP<Time> timerA = TimeMonitor::getNewCounter ("Timer A");
    RCP<Time> timerB = TimeMonitor::getNewCounter ("Timer B");
    RCP<Time> timerC = TimeMonitor::getNewCounter ("Timer C");

    // Timers A and B get nonzero elapsed times and call counts on all
    // (MPI) processes.  Timer C never gets started, so it should have
    // a zero elapsed time and zero call count on all processes.

    // The actual number of operations in the loop is proportional to
    // the cube of the loop length.  Adjust the quantities below as
    // necessary to ensure the timer reports a nonzero elapsed time
    // for each of the invocations.
    const size_t timerA_loopLength = 200;
    const size_t timerB_loopLength = 250;

    // Timer A gets a call count of 3.
    for (size_t k = 0; k < 3; ++k) {
      TimeMonitor monitorA (*timerA);
      slowLoop (size_t (timerA_loopLength));
    }
    // Timer B gets a call count of 1.
    {
      TimeMonitor monitorB (*timerB);
      slowLoop (size_t (timerB_loopLength));
    }

    const bool alwaysWriteLocal = false; // the default
    const bool writeGlobalStats = true;  // the default
    const bool writeZeroTimers = false;  // NOT the default
    TimeMonitor::summarize (oss, alwaysWriteLocal, writeGlobalStats,
                            writeZeroTimers);

    // Echo summarize() output to the FancyOStream out (which is a
    // standard unit test argument).  Output should only appear in
    // show-all-test-details mode.
    out << oss.str() << std::endl;

    // Timers A and B should be reported.
    if (Teuchos::GlobalMPISession::getRank() == 0) {
      size_t substr_i = oss.str().find ("Timer A");
      TEST_INEQUALITY(substr_i, std::string::npos);
      substr_i = oss.str().find ("Timer B");
      TEST_INEQUALITY(substr_i, std::string::npos);

      // Timer C should NOT be reported.
      substr_i = oss.str().find ("Timer C");
      TEST_EQUALITY(substr_i, std::string::npos);
    }
      
    // This sets up for the next unit test.
    TimeMonitor::clearCounters ();
  }


  //
  // Test whether the option to filter out zero timers works, for the
  // case that different (MPI) processes have different call counts.
  //
  TEUCHOS_UNIT_TEST( TimeMonitor, FILTER_ZERO_TIMERS_differentParallelCallCounts )
  {
    RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm ();
    const int myRank = comm->getRank ();

    // If MPI has not been initialized, turn the MPI communicator into
    // a "serial" communicator.  This may not be necessary when using
    // the Teuchos Unit Test Framework.
#ifdef HAVE_MPI
    {
      int mpiHasBeenInitialized = 0;
      MPI_Initialized (&mpiHasBeenInitialized);
      if (! mpiHasBeenInitialized) {
        comm = Teuchos::DefaultComm<int>::getDefaultSerialComm (null);
      }
    }
#endif // HAVE_MPI

    // Store the output of TimeMonitor::summarize() here.
    std::ostringstream oss;

    RCP<Time> timerA = TimeMonitor::getNewCounter ("Timer A");
    RCP<Time> timerB = TimeMonitor::getNewCounter ("Timer B");
    RCP<Time> timerC = TimeMonitor::getNewCounter ("Timer C");

    // Timer A gets a nonzero elapsed time and call count on Proc 0,
    // but a zero elapsed time and call count elsewhere.
    //
    // Timer B gets the same nonzero elapsed time and call count on
    // all MPI procs.
    //
    // Timer C gets a zero elapsed time and call count on all procs.
    //
    // The actual number of operations in the loop is proportional to
    // the cube of the loop length.  Adjust the quantities below as
    // necessary to ensure the timer reports a nonzero elapsed time
    // for each of the invocations.
    const size_t timerA_loopLength = 200;
    const size_t timerB_loopLength = 500;

    if (myRank == 0) {
      // Timer A gets a call count of 3 on Proc 0.
      for (int k = 0; k < 3; ++k) {
        TimeMonitor monitorA (*timerA);
        slowLoop (size_t (timerA_loopLength));
      }
    }

    // Timer B gets a call count of 1 on all procs.
    {
      TimeMonitor monitorB (*timerB);
      slowLoop (size_t (timerB_loopLength));
    }

    const bool alwaysWriteLocal = false; // the default
    const bool writeGlobalStats = true;  // the default
    const bool writeZeroTimers = false;  // NOT the default
    TimeMonitor::summarize (oss, alwaysWriteLocal, writeGlobalStats,
                            writeZeroTimers, Union);

    // Echo summarize() output to the FancyOStream out (which is a
    // standard unit test argument).  Output should only appear in
    // show-all-test-details mode.
    out << oss.str() << std::endl;

    if (Teuchos::GlobalMPISession::getRank() == 0) {
      // Timers A and B should both be reported.
      size_t substr_i = oss.str().find ("Timer A");
      TEST_INEQUALITY(substr_i, std::string::npos);
      substr_i = oss.str().find ("Timer B");
      TEST_INEQUALITY(substr_i, std::string::npos);
      
      // Timer C should NOT be reported.
      substr_i = oss.str().find ("Timer C");
      TEST_EQUALITY(substr_i, std::string::npos);
    }
      
    // This sets up for the next unit test (if there is one).
    TimeMonitor::clearCounters ();
  }

  //
  // Test option to disregard missing timers from processes in TimeMonitor::summarize().
  //
  TEUCHOS_UNIT_TEST( TimeMonitor, IgnoreMissingTimers )
  {
    RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm ();
    const int myRank = comm->getRank ();

#ifdef HAVE_MPI
    {
      int mpiHasBeenInitialized = 0;
      MPI_Initialized (&mpiHasBeenInitialized);
      if (! mpiHasBeenInitialized) {
        comm = Teuchos::DefaultComm<int>::getDefaultSerialComm (null);
      }
    }
#endif // HAVE_MPI

    // Store the output of TimeMonitor::summarize() here.
    std::ostringstream oss;

    std::string timerName="Timer Z";
    // Timer Z appears on all PIDs except 0.
    {
      switch (myRank) {
        case 1 :
          break;
        case 0 :
          {
          RCP<Time> timer = TimeMonitor::getNewCounter (timerName);
          TimeMonitor monitor (*timer);
          sleep(1);
          }
          break;
        case 2 :
          {
          RCP<Time> timer = TimeMonitor::getNewCounter (timerName);
          TimeMonitor monitor (*timer);
          sleep(1);
          }
          break;
        default :
          {
          RCP<Time> timer = TimeMonitor::getNewCounter (timerName);
          TimeMonitor monitor (*timer);
          sleep(1);
          }
    }
    }

    //////////////////////////////////////////////////////////////
    // test two versions of summarize with default behavior
    //////////////////////////////////////////////////////////////

    //version 1, comm provided
    const bool alwaysWriteLocal = false;
    const bool writeGlobalStats = true;
    const bool writeZeroTimers = false;
    bool ignoreMissingTimers = false;        // the default
    std::string filter = "";
    TimeMonitor::summarize (comm.ptr(), oss, alwaysWriteLocal, writeGlobalStats,
                            writeZeroTimers, Union, filter, ignoreMissingTimers);

    // Echo summarize() output to the FancyOStream out (which is a
    // standard unit test argument).  Output should only appear in
    // show-all-test-details mode.
    out << oss.str() << std::endl;

    if (comm->getSize() > 1) {
      // The min should be 0
      if (Teuchos::GlobalMPISession::getRank() == 0) {
        size_t substr_i = oss.str().find ("0 (0)");
        TEST_INEQUALITY(substr_i, std::string::npos);
      }
    }

    //version 2, no comm provided
    std::ostringstream oss2;
    TimeMonitor::summarize (oss2, alwaysWriteLocal, writeGlobalStats,
                            writeZeroTimers, Union, filter, ignoreMissingTimers);
    out << oss2.str() << std::endl;
    if (comm->getSize() > 1) {
      // The min should be 0
      if (Teuchos::GlobalMPISession::getRank() == 0) {
        size_t substr_i = oss2.str().find ("0 (0)");
        TEST_INEQUALITY(substr_i, std::string::npos);
      }
    }

    //////////////////////////////////////////////////////////////
    // test two versions of summarize with *non* default behavior
    //////////////////////////////////////////////////////////////
    //version 1, comm provided
    ignoreMissingTimers = true;              // NOT the default
    std::ostringstream oss3;
    TimeMonitor::summarize (comm.ptr(), oss3, alwaysWriteLocal, writeGlobalStats,
                            writeZeroTimers, Union, filter, ignoreMissingTimers);
    out << oss3.str() << std::endl;

    // The min should be different from 0
    size_t substr_i = oss3.str().find ("0 (0)");
    TEST_EQUALITY(substr_i, std::string::npos);

    //version 2, no comm provided
    std::ostringstream oss4;
    TimeMonitor::summarize (oss4, alwaysWriteLocal, writeGlobalStats,
                            writeZeroTimers, Union, filter, ignoreMissingTimers);
    out << oss4.str() << std::endl;
    // The min should be different from 0
    substr_i = oss4.str().find ("0 (0)");
    TEST_EQUALITY(substr_i, std::string::npos);

    // This sets up for the next unit test (if there is one).
    TimeMonitor::clearCounters ();
  }

  //
  // Test that Time::start() and Time::stop() call Kokkos::fence(),
  // if the option to do that is enabled.
  //
  #ifdef HAVE_TEUCHOS_TIMER_KOKKOS_FENCE
  namespace KokkosFenceCounter
  {
    static int numFences;

    void reset()
    {
      numFences = 0;
    }

    void begin_fence_callback(const char*, const uint32_t deviceId, uint64_t*) {
      using namespace Kokkos::Tools::Experimental;
      // Only count global device fences on the default space. Otherwise fences
      // could be counted multiple times, depending on how many backends are enabled.
      DeviceType fenceDevice = identifier_from_devid(deviceId).type;
      DeviceType defaultDevice = DeviceTypeTraits<Kokkos::DefaultExecutionSpace>::id;
      if(fenceDevice == defaultDevice)
        numFences++;
    }
  }

  TEUCHOS_UNIT_TEST( TimeMonitor, CheckTimerKokkosFences )
  {
    // This test doesn't care about the comm size or rank because Kokkos
    // fences and profiling is purely local to each rank.
    //
    // Set up the fence counter (reset count to 0 and set the callback)
    KokkosFenceCounter::reset();
    Kokkos::Tools::Experimental::set_begin_fence_callback(KokkosFenceCounter::begin_fence_callback);
    int fenceCountAfterStart = 0;
    int fenceCountAfterStop = 0;

    {
      RCP<Time> timer = TimeMonitor::getNewCounter ("Timer XYZ");
      TimeMonitor monitor (*timer);
      // Timer has started; record current fence count
      fenceCountAfterStart = KokkosFenceCounter::numFences;
    }
    // Timer has stopped; record again
    fenceCountAfterStop = KokkosFenceCounter::numFences;
    TEST_EQUALITY(fenceCountAfterStart, 1);
    TEST_EQUALITY(fenceCountAfterStop, 2);

    // This sets up for the next unit test (if there is one).
    TimeMonitor::clearCounters ();
  }
  #endif

} // namespace Teuchos
