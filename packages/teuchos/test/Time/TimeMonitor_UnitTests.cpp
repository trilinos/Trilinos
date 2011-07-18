/*
// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER
*/

#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_UnitTestHarness.hpp"

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


  // Slowly compute the n-th Fibonacci number.  This gives timers
  // something to time.  Be careful not to make n too large, or
  // you'll run out of stack space.
  int 
  fib (const int n) 
  {
    if (n <= 0)
      return 0;
    else if (n == 1)
      return 1;
    else
      // You should never compute the n-th Fibonacci number like
      // this.  This is exponentially slow in n.  The point of using
      // a slow algorithm like this is to exercise timers.
      return fib(n-1) + fib(n-2);
  }

  // Do a number of arithmetic operations proportional to n^3, in
  // order to have something to time.  Unlike fib() above, this
  // function shouldn't take up a lot of stack space.
  double
  slowLoop (const size_t n)
  {
    const double inc = 1 / static_cast<double> (n);
    double x = 1.0;

    for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < n; ++j)
	for (size_t k = 0; k < n; ++k)
	  x = x + inc;

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
    func_time_monitor1 ();
    std::ostringstream oss;
    TimeMonitor::summarize (oss);

    // Echo summarize() output to the FancyOStream out (which is a
    // standard unit test argument).  Output should only appear in
    // show-all-test-details mode.
    out << oss.str() << std::endl;

    const size_t substr_i = oss.str().find ("FUNC_TIME_MONITOR1");
    TEST_INEQUALITY(substr_i, std::string::npos);

    // This sets up for the next unit test.
    TimeMonitor::clearTimers ();
  }

  //
  // TimeMonitor nested timers test: create two timers on all (MPI)
  // processes, use the second inside the scope of the first, and make
  // sure that TimeMonitor::summarize() reports both timers.
  //
  TEUCHOS_UNIT_TEST( TimeMonitor, FUNC_TIME_MONITOR_tested  )
  {
    func_time_monitor2 ();
    std::ostringstream oss;
    TimeMonitor::summarize (oss);

    // Echo summarize() output to the FancyOStream out (which is a
    // standard unit test argument).  Output should only appear in
    // show-all-test-details mode.
    out << oss.str() << std::endl;

    const size_t substr_i = oss.str().find ("FUNC_TIME_MONITOR2");
    TEST_INEQUALITY(substr_i, std::string::npos);
    const size_t substr_inner_i = oss.str().find ("FUNC_TIME_MONITOR2_inner");
    TEST_INEQUALITY(substr_inner_i, std::string::npos);

    // This sets up for the next unit test.
    TimeMonitor::clearTimers ();
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
      if (! mpiHasBeenInitialized)
	comm = Teuchos::DefaultComm<int>::getDefaultSerialComm (null);
    }
#endif // HAVE_MPI

    // Store the output of TimeMonitor::summarize() here.
    std::ostringstream oss;

    // Timer A gets created on all MPI processes.
    // Timer B only gets created on Proc 0.
    RCP<Time> timerA = TimeMonitor::getNewCounter ("Timer A");
    RCP<Time> timerB;
    if (myRank == 0)
      timerB = TimeMonitor::getNewCounter ("Timer B");
    else
      timerB = null; // True anyway, but I want to make this explicit.

    // The actual number of operations in the loop is proportional to
    // the cube of the loop length.  Adjust the quantities below as
    // necessary to ensure the timer reports a nonzero elapsed time
    // for each of the invocations.
    const size_t timerA_loopLength = 150;
    const size_t timerB_loopLength = 200;

    // Timer A gets a call count of 3.
    for (size_t k = 0; k < 3; ++k)
      {
	TimeMonitor monitorA (*timerA);
	slowLoop (size_t (timerA_loopLength));
      }
    if (myRank == 0)
      { // Timer B gets a call count of 1.
	TimeMonitor monitorB (*timerB);
	slowLoop (size_t (timerB_loopLength));
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
    size_t substr_i = oss.str().find ("Timer A");
    TEST_INEQUALITY(substr_i, std::string::npos);
    if (numProcs > 1)
      {
	substr_i = oss.str().find ("Timer B");
	TEST_EQUALITY(substr_i, std::string::npos);
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
    substr_i = ossUnion.str().find ("Timer A");
    TEST_INEQUALITY(substr_i, std::string::npos);
    substr_i = ossUnion.str().find ("Timer B");
    TEST_INEQUALITY(substr_i, std::string::npos);

    // This sets up for the next unit test.
    TimeMonitor::clearTimers ();
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
    const size_t timerA_loopLength = 500;
    const size_t timerB_loopLength = 1000;

    // Timer A gets a call count of 3.
    for (size_t k = 0; k < 3; ++k)
      {
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
    size_t substr_i = oss.str().find ("Timer A");
    TEST_INEQUALITY(substr_i, std::string::npos);
    substr_i = oss.str().find ("Timer B");
    TEST_INEQUALITY(substr_i, std::string::npos);

    // Timer C should NOT be reported.
    substr_i = oss.str().find ("Timer C");
    TEST_EQUALITY(substr_i, std::string::npos);

    // This sets up for the next unit test.
    TimeMonitor::clearTimers ();
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
      if (! mpiHasBeenInitialized)
	comm = Teuchos::DefaultComm<int>::getDefaultSerialComm (null);
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

    if (myRank == 0)
      {
	// Timer A gets a call count of 3 on Proc 0.
	for (int k = 0; k < 3; ++k)
	  {
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

    // Timers A and B should both be reported.
    size_t substr_i = oss.str().find ("Timer A");
    TEST_INEQUALITY(substr_i, std::string::npos);
    substr_i = oss.str().find ("Timer B");
    TEST_INEQUALITY(substr_i, std::string::npos);

    // Timer C should NOT be reported.
    substr_i = oss.str().find ("Timer C");
    TEST_EQUALITY(substr_i, std::string::npos);

    // This sets up for the next unit test (if there is one).
    TimeMonitor::clearTimers ();
  }


} // namespace Teuchos
