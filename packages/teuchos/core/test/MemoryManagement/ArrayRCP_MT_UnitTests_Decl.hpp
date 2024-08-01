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

#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include <vector>
#include <thread>

namespace {

using Teuchos::ArrayRCP;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::null;
using Teuchos::DanglingReferenceError;

// thread method used by unit test mtArrayRCPMultipleReads below
static void read_arrayrcp_in_thread(ArrayRCP<int> shared_arrayrcp,
  int expectedValue, std::atomic<int> & countErrors) {
  // spin lock all threads until released by the main thread
  while (!ThreadTestManager::s_bAllowThreadsToRun) {}
  for( int n = 0; n < 1000; ++n) {
    // test the iterators
    for (ArrayRCP<int>::const_iterator iter = shared_arrayrcp.begin();
      iter < shared_arrayrcp.end(); ++iter) {
      // test reading a value
      int readAValue = shared_arrayrcp[0];
      // make sure the value is correct and log anything wrong
      if (readAValue != expectedValue) {
        ++countErrors;
      }
      // now check using the iterator
      int readAValueByIterator = *iter;
      // make sure the value is correct and log anything wrong
      if (readAValueByIterator != expectedValue) {
        ++countErrors;
      }
    }
  }
}

// RCP Thread Safety Unit Test: mtArrayRCPMultipleReads
//
// Purpose:
//   Sanity Check: Validate that the class is working - this was not expected
//   to have any trouble once the RCP class was made thread and no issues
//   were found.
//
// Description:
//   Creates an ArrayRCP<int>, sets all the values to an arbitrary known value,
//   then shares it to several threads which all read the ArrayRCP
//   at the same time and validate the read works. This tests both using
//   the iterators to cycle through the array and the actually reading of the
//   elements.
//
// Solution to the Problem:
//   Sanity Check
//
// Demonstration of Problem:
//   Sanity Check
TEUCHOS_UNIT_TEST( ArrayRCP, mtArrayRCPMultipleReads )
{
  const int numThreads = TEUCHOS_THREAD_SAFE_UNIT_TESTS_THREADS_USED;
  const int numTests = NUM_TESTS_TO_RUN;
  const int setValue = 67359487; // arbitrary
  const int arraySize = 10; // arbitrary
  std::atomic<int> countErrors(0); // atomic counter to log errors
  try {
    for (int testCycle = 0; testCycle < numTests; ++testCycle) {
      std::vector<std::thread> threads;
      // set up all threads to be spin locked
      ThreadTestManager::s_bAllowThreadsToRun = false;
      // create an ArrayRCP to share between the threads
      ArrayRCP<int> shared_arrayrcp(arraySize, setValue); // some array
      // Send the ArrayRCP to all the threads
      for (int i = 0; i < numThreads; ++i) {
        threads.push_back( std::thread(read_arrayrcp_in_thread,
          shared_arrayrcp, setValue, std::ref(countErrors)));
      }
      // let the threads run
      ThreadTestManager::s_bAllowThreadsToRun = true;
      // join them
      for (unsigned int i = 0; i < threads.size(); ++i) {
        threads[i].join();
      }
      convenience_log_progress(testCycle, numTests);	// this is just output
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
  TEST_EQUALITY(countErrors, 0);
}

// thread method used by unit test mtRCPofArrayRCPMultipleReads below
// this variant reads an RCP<ArrayRCP<int>> instrad of an ArrayRCP<int>
// Note the difference is the threads now directly access the same ArrayRCP<int>
// object rather than copies of it.
static void read_rcp_of_arrayrcp_in_thread(RCP<ArrayRCP<int>>
  shared_rcp_of_arrayrcp, int expectedValue, std::atomic<int> & countErrors) {
  while (!ThreadTestManager::s_bAllowThreadsToRun) {}
  for( int n = 0; n < 1000; ++n) {
    for (ArrayRCP<int>::const_iterator iter = shared_rcp_of_arrayrcp->begin();
      iter < shared_rcp_of_arrayrcp->end(); ++iter) {
      // test reading a value
      int readAValue = (*shared_rcp_of_arrayrcp)[0];
      // make sure the value is correct and log anything wrong
      if (readAValue != expectedValue) {
        ++countErrors;
      }
      // now check using the iterator
      int readAValueByIterator = *iter;
      // make sure the value is correct and log anything wrong
      if (readAValueByIterator != expectedValue) {
        ++countErrors;
      }
    }
  }
}

// RCP Thread Safety Unit Test: mtRCPofArrayRCPMultipleReads
//
// Purpose:
//   Sanity Check: Similar to prior mtArrayRCPMultipleReads test.
//
// Description:
//   Same as mtArrayRCPMultipleReads except we pass an RCP<Array<int>>
//   instead of an Array<int>
//
// Solution to the Problem:
//   Sanity Check
//
// Demonstration of Problem:
//   Sanity Check
TEUCHOS_UNIT_TEST( ArrayRCP, mtRCPofArrayRCPMultipleReads )
{
  const int numThreads = TEUCHOS_THREAD_SAFE_UNIT_TESTS_THREADS_USED;
  const int numTests = NUM_TESTS_TO_RUN;
  const int setValue = 67359487; // arbitrary
  const int arraySize = 10; // arbitrary
  std::atomic<int> countErrors(0); // atomic counter to log errors
  try {
    for (int testCycle = 0; testCycle < numTests; ++testCycle) {
      std::vector<std::thread> threads;
      // set up all threads to be spin locked
      ThreadTestManager::s_bAllowThreadsToRun = false;
      // create an RCP<ArrayRCP> to share between the threads
      RCP<ArrayRCP<int>> shared_rcp_of_arrayrcp =
        rcp(new ArrayRCP<int>(arraySize, setValue)); // some array
      // Send the RCP<ArrayRCP<int>> to all the threads
      for (int i = 0; i < numThreads; ++i) {
        threads.push_back( std::thread(read_rcp_of_arrayrcp_in_thread,
          shared_rcp_of_arrayrcp, setValue, std::ref(countErrors)) );
      }
      // let the threads run
      ThreadTestManager::s_bAllowThreadsToRun = true;
      // join them
      for (unsigned int i = 0; i < threads.size(); ++i) {
        threads[i].join();
      }
      convenience_log_progress(testCycle, numTests);	// this is just output
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
  TEST_EQUALITY(countErrors, 0);
}

} // end namespace
