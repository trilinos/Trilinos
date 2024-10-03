// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"
#ifdef HAVE_TEUCHOS_MPI
#  include "Teuchos_DefaultMpiComm.hpp"
#endif // HAVE_TEUCHOS_MPI
#include <thread>
#include <chrono>

// Check that timers are sync'd
TEUCHOS_UNIT_TEST(TimeMonitor, SyncTimeMonitor) {
  using Teuchos::SyncTimeMonitor;
  using Teuchos::Comm;
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::REDUCE_MAX;
  using Teuchos::reduceAll;
  using Clock = std::chrono::high_resolution_clock;

  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm ();

  const int numProcs = comm->getSize ();
  TEST_ASSERT( numProcs > 1 );
  if (numProcs < 4) {
    out << "This test requires at least 4 MPI processes." << std::endl;
    return;
  }

  const int myRank = comm->getRank();

  double time;
  {
    Clock::time_point start_time = Clock::now();
    {
      SyncTimeMonitor timer(*Teuchos::TimeMonitor::getNewTimer(std::string("myTimer")), comm.ptr());
      // sleep a second on rank 1 only
      if (myRank == 1)
        std::this_thread::sleep_for (std::chrono::seconds(1));
    }
    time = std::chrono::duration_cast<std::chrono::duration<double>>(Clock::now() - start_time).count();
  }

  std::ostringstream out1;
  Teuchos::TimeMonitor::summarize(out1,false,true,false,Teuchos::Union,"",true);
  std::string outStr = out1.str();
  int test = (outStr.find("myTimer") != std::string::npos);
  // Check that we took at least a second.
  int test2 = (time >= 1.0);
  test = std::min(test, test2);

  int gblTest = false; // output argument
  reduceAll<int, int> (*comm, REDUCE_MAX, test, outArg (gblTest));
  TEST_EQUALITY(gblTest, 1);

}

// Check that sync'd timers do not hang execution.
TEUCHOS_UNIT_TEST(TimeMonitor, HangingSyncTimeMonitor) {
  using Teuchos::SyncTimeMonitor;
  using Teuchos::Comm;
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::REDUCE_MAX;
  using Teuchos::reduceAll;

  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm ();

  const int numProcs = comm->getSize ();
  TEST_ASSERT( numProcs > 1 );
  if (numProcs < 4) {
    out << "This test requires at least 4 MPI processes." << std::endl;
    return;
  }

  const int myRank = comm->getRank();
  int test = false;

  try {
    // Setup up a sync'd timer and raise an exception on one rank.
    {
      SyncTimeMonitor timer(*Teuchos::TimeMonitor::getNewTimer(std::string("myTimer")), comm.ptr());
      if (myRank == 1)
        throw std::runtime_error("Test");
    }
    test = true;
  } catch (const std::runtime_error& e) {
    test = (myRank == 1);
  }

  int gblTest = false; // output argument
  reduceAll<int, int> (*comm, REDUCE_MAX, test, outArg (gblTest));
  TEST_EQUALITY(gblTest, 1);

}
