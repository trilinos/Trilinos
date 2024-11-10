// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// These unit tests are used for both a Nightly version and a Basic version

#include "General_MT_UnitTests.hpp"

#include "Teuchos_ArrayView.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include <vector>
#include <thread>

namespace {

using Teuchos::ArrayView;
using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Array;
using Teuchos::ArrayRCP;
using Teuchos::RangeError;
using Teuchos::DanglingReferenceError;

// utlity method used by unit test mtArrayViewMultipleReads below
// check the iterators don't do anything bad
static void read_arrayview_in_thread(RCP<ArrayView<int>> shared_arrayview,
  int expectedValue, std::atomic<int> & countErrors) {
  // spin lock the threads
  while (!ThreadTestManager::s_bAllowThreadsToRun) {}
  for( int n = 0; n < 1000; ++n) {
    for (ArrayView<int>::iterator iter = shared_arrayview->begin();
      iter < shared_arrayview->end(); ++iter) {
      int readAValueByIterator = *iter;
      // make sure the value is correct and log anything wrong
      if (readAValueByIterator != expectedValue) {
        ++countErrors;
      }
    }
  }
}

// RCP Thread Safety Unit Test: mtArrayViewMultipleReads
//
// Purpose:
//   Sanity Check: Validate that the class is working - this was not expected
//   to have any trouble once the RCP class was made thread and no issues
//   were found.
//
// Description:
//   Creates an Array<int>, sets all the values to an arbitrary known value,
//   then create an RCP<ArrayView<int>> of that array,
//   then share it to several threads which all read the ArrayView
//   at the same time and validate the read works. This tests both using
//   the iterators to cycle through the array and the actually reading of the
//   elements. This mirrors the Array test (which will fail without mutex
//   protection) but ArrayView is ok because the ArrayView begin does not have
//   any mutable behavior (it is true const).
//
// Solution to the Problem:
//   Sanity Check
//
// Demonstration of Problem:
//   Sanity Check
TEUCHOS_UNIT_TEST( ArrayView, mtArrayViewMultipleReads )
{
  const int numThreads = TEUCHOS_THREAD_SAFE_UNIT_TESTS_THREADS_USED;
  const int numTests = NUM_TESTS_TO_RUN;
  const int setValue = 67359487; // arbitrary
  const int arraySize = 10; // arbitrary
  std::atomic<int> countErrors(0); // atomic counter to log errors
  for (int testCycle = 0; testCycle < numTests; ++testCycle) {
    try {
      std::vector<std::thread> threads;
      ThreadTestManager::s_bAllowThreadsToRun = false;
      Array<int> array(arraySize, setValue); // some array
      RCP<ArrayView<int>> arrayview_rcp = rcp(new ArrayView<int>(array));

      for (int i = 0; i < numThreads; ++i) {
        threads.push_back( std::thread(read_arrayview_in_thread,
          arrayview_rcp, setValue, std::ref(countErrors)) );
      }

      ThreadTestManager::s_bAllowThreadsToRun = true;   // let the threads run
      for (unsigned int i = 0; i < threads.size(); ++i) {
        threads[i].join();
      }
    }
    TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);

    convenience_log_progress(testCycle, numTests);	// this is just output
  }

  // right now this test is just looking for trouble - so we don't actually
  // have a verification - a failed test would be corrupted memory for example
  TEST_EQUALITY_CONST(0, 0);
}

// this test is only meaningful in DEBUG and would crash in RELEASE
// with undefined behaviors
#ifdef TEUCHOS_DEBUG

// this utility method is defined to create and delete memory
// the purpose of this thread is just to put pressued on memory
// so we can study whether the classes are thread safe
// Note this test closely mirrors the Ptr unit test which detects scrambled
// memory events.
static void scramble_memory(int scrambleValue, int testArraySize,
  int finishWhenThisThreadCountCompletes) {
  while (!ThreadTestManager::s_bAllowThreadsToRun) {}
  // the idea here is to try and fill any heap holes with new int allocations
  while (true) {
    // hard coded this as for thread debugging I didn't want to have extra array
    // methods running while investigating the main operations
    #define ARRAY_SCRAMBLE_SIZE 100
    std::vector<int> * tempPtrArray[ARRAY_SCRAMBLE_SIZE];
    for (int n = 0; n < ARRAY_SCRAMBLE_SIZE; ++n) {
      // if the scramble thread does not allocate std::vector chunks identical
      // to the main thread it won't have any chance to trigger the scrambled
      // events.
      tempPtrArray[n] = new std::vector<int>(testArraySize, scrambleValue);
    }
    for (int n = 0; n < ARRAY_SCRAMBLE_SIZE; ++n) {
      delete tempPtrArray[n];
    }
    if (ThreadTestManager::s_countCompletedThreads >=
      finishWhenThisThreadCountCompletes) {
      break;
    }
  }
}

// note this test mirrors the Ptr test - the mechanisms we are considering
// here are essentially equivalent. When a weak RCP is raced, it can indicate
// that it is not dangling, and then be read on junk memory because the delete
// happens right after the dangling check
static void share_arrayview_to_threads(ArrayView<int> shared_arrayview,
  int theTestValue, Cycle_Index_Tracker & index_tracker) {
  while (!ThreadTestManager::s_bAllowThreadsToRun) {}
  int cycle = 0;
  try {
    while (true) {
      bool bCheckStatus = ThreadTestManager::s_bMainThreadSetToNull;
      // this may return junk data if the new heap allocation jumped in
      // any of the member values could be junk
      int tryToReadAValue = shared_arrayview[0];
      index_tracker.trackCycle = cycle;
      if (tryToReadAValue != theTestValue) {
        // if we get here we had an ok from the dangling reference check, but
        // then memory was deleted and reallocated to a new value by the
        // scramble thread - a rare but possible condition
        index_tracker.scambledMemory = cycle;
      }

      if (bCheckStatus) {
        index_tracker.unknownError = cycle;
        // when bCheckStatus true it means we started the loop after the main
        // rcp was set null - we should have thrown a DanglingReference by now
        break;
      }
      ++cycle;
    }
  }
  catch (DanglingReferenceError&) {
    // loop ends - we got the dangling reference
    index_tracker.danglingReference = cycle;
  }
  catch (...) {
    std::cout << std::endl << "Unknown and unhandled exception!" << std::endl;
  }

  ++ThreadTestManager::s_countCompletedThreads;
}

// RCP Thread Safety Unit Test: mtArrayViewDangling
//
// Purpose:
//   To understand this test is may be worth first looking at
//   the Ptr unit test mtPtrDangling which studies a similar mechanism.
//   This test is partly a sanity check to make sure dangling references are
//   properly handled.
//
// Description:
//   In this test an ArrayView is created for an ArrayRCP, then passed
//   to several threads. The main thread kills the ArrayRCP and the subthreads
//   all get left with dangling weak references. The test collects data on
//   whether the subthreads properly process this.
//   There is one outstanding issue here where a weak RCP can think it is valid
//   but then read invalid memory. So this means Debug can detect problems
//   most of the time but may sometimes do an undefined behavior. Currently
//   there is no fix in place. The mtPtrDangling successfully demonstrates
//   these bad memory reads and detects them. This should happen here as well
//   but for some reason they never seem to occur - maybe due to the larger
//   size of the memory structures making it hard to have some other thread
//   jump in an replace the values. So this test shows that dangling references
//   are detected but does not show the full picture the way mtPtrDangling does.
//
// Solution to the Problem:
//   The dangling references work as expected due to the other RCP changes.
//
// Demonstration of Problem:
//   To see the weak RCP accessing bad memory without knowing it, currently
//   use the mtPtrDangling test as this one doesn't seem to trigger it and it's
//   not clear why.
TEUCHOS_UNIT_TEST( ArrayView, mtArrayViewDangling )
{
  const int numThreads = TEUCHOS_THREAD_SAFE_UNIT_TESTS_THREADS_USED;
  const int numTests = NUM_TESTS_TO_RUN;
  const int theTestValue = 66387; // some value
  const int scrambleValue = 572778; // some other value
  const int testArraySize = 3;

  // we want to count when it's not trivial (first cycle or last cycle)
  int countDanglingReferences = 0;
  int scrambledMemoryEvents = 0;
  int unknownErrors = 0;
  // 0 is the scrambling thread doing constant new/delete.
  // The rest are the reader threads looking for troubles
  int finishWhenThisThreadCountCompletes = numThreads - 1;
  for (int testCycle = 0; testCycle < numTests; ++testCycle) {
    try {
      ThreadTestManager::s_countCompletedThreads = 0;

      // first make an arrayrcp which we will kill later
      ArrayRCP<int> arrayrcp = arcp(rcp(
        new std::vector<int>(testArraySize, theTestValue)));
        // now make an ArrayView which has a reference to the arrayrcp
      ArrayView<int> shared_arrayview = arrayrcp();
      // threads will start spin locked
      ThreadTestManager::s_bAllowThreadsToRun = false;
      // this will be used to tell the threads when we have killed the RCP
      ThreadTestManager::s_bMainThreadSetToNull = false;
      // used to track errors in the sub threads
      Cycle_Index_Tracker index_tracker[numThreads];
      // create the threads
      std::vector<std::thread> threads;
      for (int i = 0; i < numThreads; ++i) {
        switch(i) {
          case 0:
          {
            // the first thread is special and just puts pressure on memory
            // but allocating and deleting - tried some different combinations
            // here but could never get this thread to jump in to the released
            // memory spot of the ArrayRCP as in mtPtrDangling.
            // The larger data structure probably is making this harder to
            // demonstrate.
            threads.push_back(std::thread(scramble_memory, scrambleValue,
              testArraySize, finishWhenThisThreadCountCompletes));
          }
          break;
          default:
          {
            // These threads all just read the ArrayView and will process
            // dangling reference exceptions when the ArrayRCP is killed by this
            // main thread.
            threads.push_back(std::thread(share_arrayview_to_threads,
              shared_arrayview, theTestValue, std::ref(index_tracker[i])));
          }
          break;
        }
      }
      // let the threads start running
      ThreadTestManager::s_bAllowThreadsToRun = true;
      // spin lock until we have confirmed the sub threads did something
      while (index_tracker[1].trackCycle < 1) {}
      // the ArrayRCP becomes invalid and the ArrayView types all lose their
      // valid object - now we start getting dangling references
      arrayrcp = null;
      ThreadTestManager::s_bMainThreadSetToNull = true;   // tell the threads
      // join the threads
      for (unsigned int i = 0; i < threads.size(); ++i) {
        threads[i].join();
      }
      // collect all the errors
      for (unsigned int i = 0; i < threads.size(); ++i) {
        if (index_tracker[i].danglingReference != -1) {
          ++countDanglingReferences; // this should always happen
        }
        if (index_tracker[i].scambledMemory != -1 ) {
          // this was expected but does not occur
          // in the mtPtrDangling these events are detected
          // presently it's not clear why this test could not also demonstrate
          // this shortcoming of weak RCP in the present setup.
          ++scrambledMemoryEvents;
        }
        if (index_tracker[i].unknownError != -1 ) {
          ++unknownErrors; // this is not expected and never occurs
        }
      }
    }
    TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
    convenience_log_progress(testCycle, numTests);	// this is just output
  }

  // verify we got all the dangling references
  // except for thread 0 (scramble thread) we should be getting 1 exception
  // for every run in each thread.
  int requiredDanglingReferenceCount = (numThreads-1) * numTests;
  bool bDanglingReferenceDetectionCountIsOK = (countDanglingReferences ==
    requiredDanglingReferenceCount);

  // if the dangling exception count was off log some information
  if( !bDanglingReferenceDetectionCountIsOK ) {
    std::cout << std::endl << "Detected " << countDanglingReferences <<
      " Dangling References but should have found " <<
      requiredDanglingReferenceCount << "." << std::endl;
  }
  else {
    // if everything is ok log the info along with scrambled memory events
    // scrambles is always 0 here but I would not be surprised if it could
    // be non zero. Currently it is not a factor for whether the test will fail.
    std::cout << "Danglers: " << countDanglingReferences << " Scrambles: " <<
      scrambledMemoryEvents << " ";
  }

  // this is not expected to occur and was not ever observed
  if (unknownErrors != 0) {
    std::cout << std::endl << "Detected " << unknownErrors <<
      " dangling references were missed which should have been detected."
      << std::endl;
  }
  // pass or fail the test
  TEST_ASSERT( bDanglingReferenceDetectionCountIsOK )
  // if this ever hits it would be unexpected
  TEST_EQUALITY_CONST(unknownErrors, 0);
}

#endif // TEUCHOS_DEBUG

} // end namespace
