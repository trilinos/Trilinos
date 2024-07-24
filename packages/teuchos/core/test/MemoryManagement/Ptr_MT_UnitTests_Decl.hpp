// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// These unit tests are used for both a Nightly version and a Basic version

// this test is only meaningful in DEBUG and would crash in RELEASE
// with undefined behavior. This is because the test involves debug checks
// to detect weak ptrs and in release these errors are ignored. So here we
// are checking whether debug code can safely detectly badly written code.
#include "Teuchos_ConfigDefs.hpp" // get TEUCHOS_DEBUG

#ifdef TEUCHOS_DEBUG

#include "General_MT_UnitTests.hpp"
#include "Teuchos_Ptr.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include <vector>
#include <thread>
#include <atomic>

namespace {

using Teuchos::Ptr;
using Teuchos::RCP;
using Teuchos::DanglingReferenceError;
using Teuchos::null;
using Teuchos::rcp;
using Teuchos::ptrFromRef;
using Teuchos::rcpFromPtr;

// method used by unit test mtPtrDangling below.
// the thread reads the shared Ptr<int> which has been release by the
// main thread. The weak RCP is intended to detect when this is read after
// being released, which is a programmer error.
// The thread also puts pressue on memory by allocating/deleting ints.
static void share_ptr_to_threads(Ptr<int> shared_ptr, int theTestValue,
  Cycle_Index_Tracker & index_tracker) {
  // spin lock the threads until release by the main thread
  while (!ThreadTestManager::s_bAllowThreadsToRun) {}
  int cycle = 0;
  try {
    // If there is lots of competition for threads setting this to some
    // safety limit of counts may fail because another thread was held up.
    // So looping while(true) may be the cleanest and then we just
    // time out if something goes wrong.
    while(true) {
      // check if the main thread has released the RCP which we point to.
      bool bCheckStatus = ThreadTestManager::s_bMainThreadSetToNull;
      // keep track of which cycle we are on
      index_tracker.trackCycle = cycle;

      // Now read the ptr - there are 4 possible outcomes:
      // (1) the ptr debug check returns dangling and a proper throw is
      //     detected - in this case we are certain of our result
      // (2) the ptr debug check returns valid and we can read the data
      //     (because we are lucky and the data remains valid while we use it)
      // (3) the ptr debug check returns valid, gets deleted by another
      //     thread immediately after, but we read the deleted data without
      //      knowing because it still contains the proper memory
      // (4) the ptr debug check returns valid, gets deleted by another
      //     thread immediately after, is overwriteen by another heap
      //     allocation, and we read the scrambled data without knowing
      if (*shared_ptr != theTestValue) {
        index_tracker.scambledMemory = cycle; // record the cycle of the error
      }

      // the scrambler int is trying to jump into the released memory spot
      // through a heap allocation and disrupt the ptr value
      int * pScramblerInt = new int;
      *pScramblerInt = 0; // we hope to set the dangling memory space here
      delete pScramblerInt;

      // if the main thread had released the memory before we read the ptr
      // then we should have thrown by now. So something else has gone wrong
      // and we record an unknown error (this currently does not every happen).
      if (bCheckStatus) {
        index_tracker.unknownError = cycle;
        break;
      }
      ++cycle;
    }
  }
  catch(DanglingReferenceError&) {
    // we got the dangling error as expected
    index_tracker.danglingReference = cycle;
  }
}

// RCP Thread Safety Unit Test: mtPtrDangling
//
// Purpose:
//   Validate the RCP Ptr mechanism are all thread safe.
//   Currently Ptr can detect a dangling reference if the original RCP was
//   released in another thread. However because it is based on a weak RCP
//   mechanism it can be fooled and think the Ptr was valid but then before
//   actually reading the memory, lose that value.
//   So this test will show the danlging references are processed properly
//   which happens almost every time. However occasionally the scrambled memory
//   event will occur (currently no fix) and this test is designed to detect
//   that it took place. At some point if we decide to fix this we can use this
//   test to validate it's all working.
//
// Description:
//   An RCP<int> is created and then a Ptr<int> is creaeted from the RCP which
//   maintains a weak reference to the original RCP<int>. The subthreads read
//   the Ptr<int> while the main thread releases the memory. The subthreads then
//   detect they have a dangling reference and throw an exception.
//
// Solution to the Problem:
//   There is currently no implemented solution for the debug situation where
//   an RCP class attempts to dereference a weak ptr
//
// Demonstration of Problem:
//   Running this test will provide output showing dangling reference being
//   properly detected and a few srambled memory events occuring which currently
//   does not have a fix in place.
TEUCHOS_UNIT_TEST( Ptr, mtPtrDangling )
{
  const int numThreads = TEUCHOS_THREAD_SAFE_UNIT_TESTS_THREADS_USED;
  const int numTests = NUM_TESTS_TO_RUN;
  const int theTestValue = 1454083084; // see Ptr to arbitrary value
  // we want to count when it's not trivial (first cycle or last cycle)
  int countDanglingReferences = 0; // detect attempt to access deleted RCP
  int scrambledMemoryEvents = 0; // detect scambled memory
  int unknownErrors = 0; // detect unknown errors - currently shouldn't happen
  for (int testCycle = 0; testCycle < numTests; ++testCycle) {
    try {
      // create a new int - RCP will own this int and manage its memory
      int * pInt = new int;
      // set the int to a test value - we will check for this value in threads
      *pInt = theTestValue;
      // first make an RCP
      RCP<int> shared_rcp = rcp(pInt);
      // now make a Ptr which remembers a weak reference to that RCP
      Ptr<int> shared_ptr = shared_rcp.ptr();
      // threads will start spin locked
      ThreadTestManager::s_bAllowThreadsToRun = false;
      // we have not yet deleted the RCP in this thread
      ThreadTestManager::s_bMainThreadSetToNull = false;
      // manager to keep track of events
      Cycle_Index_Tracker index_tracker[numThreads];
      // Now create the threads
      std::vector<std::thread> threads;
      for (int i = 0; i < numThreads; ++i) {
        threads.push_back(std::thread(share_ptr_to_threads, shared_ptr,
          theTestValue, std::ref(index_tracker[i])));
      }
      // let the threads run
      ThreadTestManager::s_bAllowThreadsToRun = true;
      // spin lock the main thread until the sub threads get started
      while( index_tracker[0].trackCycle < 1 ) {}
      // Now set the RCP null
      // the RCP becomes invalid and the Ptr types all lose their valid object
      shared_rcp = null;
      // tell the threads the RCP is now dead
      // This lets the threads know they 'must' detect errors on next loop
      ThreadTestManager::s_bMainThreadSetToNull = true;
      // Join all threads to completion
      for (unsigned int i = 0; i < threads.size(); ++i) {
        threads[i].join();
      }
      // count up all the errors
      for (unsigned int i = 0; i < threads.size(); ++i) {
        if (index_tracker[i].danglingReference != -1) {
          ++countDanglingReferences; // common event
        }
        if (index_tracker[i].scambledMemory != -1 ) {
          ++scrambledMemoryEvents; // happens but rarely
        }
        if (index_tracker[i].unknownError != -1 ) {
          ++unknownErrors; // not presently ever an issue
        }
      }
    }
    TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
    convenience_log_progress(testCycle, numTests);// this is just output
  }

  // verify we caught a dangler everytime
  int expectedDanglingReferences = numThreads * numTests;
  if( countDanglingReferences != expectedDanglingReferences) {
    std::cout << std::endl << "Test FAILED because only " <<
      countDanglingReferences <<
      " dangling references were detected but expected "
      << expectedDanglingReferences << "." << std::endl;
  }
  else {
    // we got the expected number of danglers so this log is the output
    // when the test succeeeds. We also log the number of scambled events.
    // At some point we may implement a fix so that this missed event does not
    // occur but that is not currently the case. The weak RCP can be tricked,
    // think it's valid, and read the memory which subsequently was deleted.
    // If we do implement a fix in the future, we can check here as scrambled
    // events should then go to 0. We would then always detect invalid memory
    // in debug mode.
    std::cout << "Danglers: " << countDanglingReferences << " Scrambles: "
      << scrambledMemoryEvents << " ";
  }

  // this is not currently an issue - it was a safety check in case something
  // unexpected ever happened in the thread loop
  if (unknownErrors != 0) {
    std::cout << std::endl << "Detected " << unknownErrors <<
      " dangling references were missed which should have been detected."
        << std::endl;
  }
  // verify we detected the expectedDanglingReferences
  TEST_ASSERT(countDanglingReferences == expectedDanglingReferences)
  // not presently an issue - this is searching for the possibility of a
  // dangling reference missed when it should have been recorded
  TEST_EQUALITY_CONST(unknownErrors, 0);
}

} // end namespace

#endif // TEUCHOS_DEBUG
