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

#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include <vector>
#include <thread>

namespace {

using Teuchos::Array;
using Teuchos::RCP;
using Teuchos::DanglingReferenceError;
using Teuchos::RangeError;

// this thread method repeatedly calls the iterator loop and reads values
// arrayType will be Array or const Array
// iteratorType will be corresponding ::iterator or ::const_iterator
template<typename arrayType, typename iteratorType>
static void thread_reads_array(arrayType shared_array, int setValue) {
  while (!ThreadTestManager::s_bAllowThreadsToRun) {}
  for( int n = 0; n < 1000; ++n) {
    for (iteratorType iter = shared_array->begin();
      iter < shared_array->end(); ++iter) {
      int readValue = *iter; // read the value
      if(readValue != setValue) {
        throw std::logic_error("Test failed to read proper array value.");
      }
    }
  }
}

// Shared method for Array and const Array unit tests
// See notes below for individual unit tests
template<typename arrayType, typename iteratorType>
void runReadArrayTest()
{
  const int numThreads = TEUCHOS_THREAD_SAFE_UNIT_TESTS_THREADS_USED;
  const int numTests = NUM_TESTS_TO_RUN;
  const int setArrayLength = 10; // somewhat arbitrary
  const int setArrayValue = 3; // arbitrary
  for (int testCycle = 0; testCycle < numTests; ++testCycle) {
    std::vector<std::thread> threads;
    // set up threads to be spin locked
    ThreadTestManager::s_bAllowThreadsToRun = false;
    // makes an array of length setArrayLength with each
    // element set to setArrayValue
    arrayType array_rcp =
      rcp(new Array<int>( setArrayLength, setArrayValue ));
    // create multiple threads which read the array
    for (int i = 0; i < numThreads; ++i) {
      threads.push_back(std::thread(
        thread_reads_array<arrayType, iteratorType>,
          array_rcp, setArrayValue));
    }
    // let the threads run
    ThreadTestManager::s_bAllowThreadsToRun = true;
    // join the threads
    for (unsigned int i = 0; i < threads.size(); ++i) {
      threads[i].join();
    }
    convenience_log_progress(testCycle, numTests);  // just output
  }
}

// RCP Thread Safety Unit Test: mtArrayMultipleReads_NonConst
//
// Purpose:
//   Demonstrates the Array class needs thread safety protection for debug mode
//   since the debug ArrayRCP objects are allocated in the begin() function.
//
// Description:
//   An Array is shared to multiple threads - they simultaneously use the
//   non-const iterator - so this calls begin(), end(), and ++iter on the
//   thread which was sufficient to demonstrate the original problem.
//   Later added reads of the iterator for more complete testing.
//   The point of this test was to validate that multiple threads can safely
//   read a shared Array. We expected it to fail originally because the begin()
//   call in Debug will set extern_arcp_ so the first strategy was to make a
//   race condition on that allocation to demonstrate we could see this problem.
//   Note that begin() is a const but the internal extern_arcp_ object is
//   mutable - that object could be changed by multiple threads simultaneously.
//
// Solution to the Problem:
//   A single debug-only mutex protects all changes to the member vec_ as well
//   as the debug ArrayRCP objects to avoid any simultaneous read/write
//   behavior.
//
// Demonstration of Original Problem:
//   Add the following define at top of this file:
//     #define REMOVE_THREAD_PROTECTION_FOR_ARRAY
//   That will remove the mutex protecting the debug mode.
TEUCHOS_UNIT_TEST( Array, mtArrayMultipleReads_NonConst )
{
  try {
    runReadArrayTest<RCP<Array<int>>, Array<int>::iterator>();
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
}

// RCP Thread Safety Unit Test: mtArrayMultipleReads_Const
//
// Purpose:
//   Similar to mtArrayMultipleReads_NonConst.
//
// Description:
//   Similar to mtArrayMultipleReads_NonConst.
//
// Solution to the Problem:
//   Similar to mtArrayMultipleReads_NonConst though the original Array form had
//   the const begin() call the non-const begin(), so this was refactored to
//   facilitate the mutex protection
//
// Demonstration of Original Problem:
//   Similar to mtArrayMultipleReads_NonConst.
TEUCHOS_UNIT_TEST( Array, mtArrayMultipleReads_Const )
{
  try {
    runReadArrayTest<const RCP<Array<int>>, Array<int>::const_iterator>();
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
}

// These tests are only meaningful in DEBUG
// They would crash with undefined behaviors in RELEASE
// All of this Array work was designed to make the debug tracking for Array
// work in a thread safe way but in release those protections don't exist.
#ifdef TEUCHOS_DEBUG

// This method runs continuously inserting setValue on shared_array and then
// removing those values. The shared_array never gets larger than maxArraySize
// The process completes when finishWhenThisThreadCountCompletes determines
// that all the reading threads have completed their operations.
// For the unit tests there will only be 1 thread calling this method.
static void call_inserts_on_array(RCP<Array<int>> shared_array, int setValue,
  int maxArraySize, int finishWhenThisThreadCountCompletes) {
  while (!ThreadTestManager::s_bAllowThreadsToRun) {}
  while (true) { // runs until read threads register completion
    const int insertCount = maxArraySize - shared_array->size();
    // insert some values
    for (int n = 0; n < insertCount; ++n) {
      shared_array->push_back(setValue); // insert ints with value setValue
    }
    // erase some values
    for (int n = 0; n < insertCount; ++n) {
      shared_array->pop_back(); // remove values so it doesn't get too big
    }
    if (ThreadTestManager::s_countCompletedThreads >=
      finishWhenThisThreadCountCompletes) {
      break;
    }
    // track cycles - used to make sure this thread is up and running before
    // we begin testing the read threads
    ++ThreadTestManager::s_countWritingThreadCycles;
  }
}

// This method runs continuously deleting and allocated memory.
// When Weak RCP checks for a valid ptr and another strong RCP deallocates,
// the memory may be lost after the check but before the actual read.
// This scramble thread can claim that memory, and replace the deallocated
// memory spot with a new scramble value. This allows the unit test to detect
// that this bad memory event has occurred. Note that currently we do not
// have a fix for this Debug Weak RCP check
static void scramble_memory(int scrambleValue,
  int finishWhenThisThreadCountCompletes) {
  while (!ThreadTestManager::s_bAllowThreadsToRun) {}
  while (true) {
    // here I used a strict C array to avoid having this code interfering with
    // any debugging of the actual Array class
    // Allocates 100 int ptrs, sets them, then deletes the memory, then repeat
    // Better than just new-delete repeated because we don't want to just
    // attack the same memory slot over and over.
    #define ARRAY_SCRAMBLE_SIZE 100
    int * tempPtrArray[ARRAY_SCRAMBLE_SIZE];
    for (int n = 0; n < ARRAY_SCRAMBLE_SIZE; ++n) {
      int * pInt = new int;
      *pInt = scrambleValue;
      tempPtrArray[n] = pInt;
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

// defined three modes of array testing
// during these tests one thread will be inserting/removing values and
// another thread will be allocating/deleting memory.
// All other threads will perform the operation defined here.
enum ArrayTest_Style {
  // Will try to read array - this currently can lead to out of range errors
  // no iteration occurs in thie mode.
  ArrayTest_DoReadOperations,

  // Iterator will trigger dangling references but this simple form may
  // not trigger dangling the first cycle just due to random luck of the draw.
  ArrayTest_TriggerDanglingWithIteration,

  // Modified test to hit dangling reference immediately on the first cycle.
  // It's a sanity check and the test happens in a deterministic way.
  ArrayTest_TriggerDanglingWithIterationFirstCycle
};

// This method will be called on 2 threads and run an operation
// defined by arrayTestStyle. While those 2 threads run there will also be
// 2 special threads running - one is inserting/removing values while another
// is constantly allocating/deleting memory.
template<class arrayType>
static void do_read_operations_on_array(RCP<arrayType> shared_array,
  int setValue, int scrambleValue, Cycle_Index_Tracker & index_tracker,
  int maxArraySize, ArrayTest_Style arrayTestStyle) {
  // spin lock until the threads are released
  while (!ThreadTestManager::s_bAllowThreadsToRun) {}
  int cycle = 0; // track the loop cycle so we can save it
  try {
    // once we catch the event this loop ends so the count should not impact
    // the time for the test. 10000 is just to have a safety exit and then the
    // test can fail due to not detecting the events.
    for (cycle = 0; cycle < 10000; ++cycle) {
      switch (arrayTestStyle) {
        case ArrayTest_DoReadOperations:
        {
          // read some values and check for scrambles. Two things can happen,
          // either we get a range error and get a valid exception or perhaps
          // we just straight out read bad memory in which case we won't get
          // an exception but will try to trap that event here and record it.
          int readValue = shared_array->at(0);
          if( readValue != setValue ) {
            // was using this to see if things were working properly - but we
            // don't necessarily expect the scrambled int to always measure
            // here - something else could be going on.
            if( readValue == scrambleValue) {
              // we detected a bad memory read and the int was set by the
              // scrambling thread. That confirms we really did read memory
              // which was allocated by another thread
              index_tracker.scambledMemory = cycle;
            }
            else {
              // we detected bad memory but it was not set by the scrambling
              // thread. Anything else could be going on to have changed this.
              // For now just do the same as above and mark we detected the
              // error. We are not recording the distinction between this and
              // the above case right now. I have left both entries for debug
              // testing.
              index_tracker.scambledMemory = cycle;
            }
          }
        }
        break;
        case ArrayTest_TriggerDanglingWithIteration:
        {
          // this is perhaps the most 'natural' way to test the danglers
          // detection using iterators. However note that this is not
          // guaranteed to trigger a dangling reference on a particular cycle.
          // shared_array()->begin() and shared_array->end() do not depend on
          // the iter arrayRCP so they won't necessarily pick up the dangling
          // reference immediately. It's possible for shared_array to be
          // changing in a state so that iter++ never gets called after it
          // gets set to shared_array->begin(). For example
          // shared_array()->end() may change to be much smaller than begin(),
          // or much larger(). Then the loop ends and we go to the next cycle.
          // That is why this test loops in the first place.
          // Another test bGuaranteeDanglingDetectionOnFirstCycle true is
          // designed to always hit dangling on the first cycle for clarity
          // and to provide a deterministic test mode.
          for (Array<int>::const_iterator iter = shared_array->begin();
            iter < shared_array->end(); ++iter) {
            // empty loop because we are testing the iterators and how they
            // behave when another thread changes the array in the middle of
            // iterating.
          }
        }
        break;
        case ArrayTest_TriggerDanglingWithIterationFirstCycle:
        {
          // this deterministic mode is designed to hit the dangling reference
          // on the very first cycle. If we got to the second cycle we didn't
          // succeed and need to report an error.
          if (cycle != 0) {
            // log the event so the test can fail
            index_tracker.missedDanglingOnFirstCycle = cycle;

            // inform other threads now that this thread is completed
            // so they can quit properly
            ++ThreadTestManager::s_countCompletedThreads;

            // this will exit the thread loop immediately
            return;
          }

          // This is the logic to make this threead always hit the dangling
          // reference on the first cycle. First we allocate an iterator.
          Array<int>::const_iterator iter = shared_array->begin();

          // determine what thread cycle the read/write thread is currently on
          int getReadWriteThreadCycles =
            ThreadTestManager::s_countWritingThreadCycles;

          // spin lock this thread until the read/write thread does at
          // least one full cycle beyond where it was 'after' the moment
          // we created the iterator.
          while (ThreadTestManager::s_countWritingThreadCycles <
            getReadWriteThreadCycles + 2) {  // empty loop
          }

          // Now call ++iter once - this MUST fail and detect the dangler
          // In this condition iter is always a weak reference to a node with
          // an invalid ptr.
          ++iter; // the test must throw here!
        }
        break;
      }
    }
  }
  catch (DanglingReferenceError&) {
    // If test throws a dangling reference error, record the cycle it occurred.
    index_tracker.danglingReference = cycle;
  }
  catch (RangeError&) {
    // If tests throws a range error, record the cycle it occurred.
    index_tracker.outOfRangeError = cycle;
  }
  catch (std::out_of_range&) {
    // We could also get std::out_of_range
    // Note that here were are counting Trilinos RangeError and STL
    // std::out_of_range as all the same - they both happen since the bad
    // memory read can be triggered at any time and it depends on exactly where
    // we were in the call stack when things happened.
    index_tracker.outOfRangeError = cycle;
  }

  // advertise we are done. When all these threads complete,
  // the push/pop thread and scambler thread will quit as well.
  ++ThreadTestManager::s_countCompletedThreads;
}

// This method will set up:
// One thread which is constantly inserting/deleting values form shared_array.
// Another thread will constantly allocate/delete
// Remaining threads will all run do_read_operations_on_array (above)
// with a read operation defined by arrayTestStyle.
// This method is used for all the unit tests except for the two simple
// read tests defined above.
bool runArrayDanglingReferenceTest( bool bUseConstVersion,
  ArrayTest_Style arrayTestStyle ) {
  // Note that 0 is the insert thread, and 1 is the scrambler thread so set
  // this to be 4 or more - 3 is not ideal as only have one thread reading.
  // Using 4 because NUM_TOTAL_CORES_USED - 4
  const int numThreads = TEUCHOS_THREAD_SAFE_UNIT_TESTS_THREADS_USED;

  // how many times to repeat the entire test
  const int numTests = NUM_TESTS_TO_RUN;

  // arbitrary - set up an array filled with this value
  const int setValue = 1;

  // arbitrary - but not the same as setValue
  // the scramble thread tries to allocate memory during the race condition
  // which occurs for weak RCP when it checks for valid memory and then loses
  // that memory before actually reading it. By setting a known value the unit
  // test can confirm it is registering that event. Exactly how we will
  // resolve that is a pending issue.
  const int scrambleValue = 12345;

  // The insert/delete thread will build the array up to this size and
  // then start to delete elements
  const int maxArraySize = 100;

  // the maximum number of errors we can detect.
  // One per thread which is doing read operations.
  // Does not count the push/pop thread or the scramble thread.
  int countTotalTestRuns = 0;

   // the actual number of dangling references we detect
  int countDetectedDanglingReferences = 0;

  // if using the test designed to hit first cycle, this tracks failures
  int countMissedFirstCycleDanglers = 0;

  // this counts the times the dangling reference missed (thought memory was
  // fine) but it was actually not and it was found to be overwritten to the
  // scrambleValue
  int countScrambledMemoryEvents = 0;

  // this counts out of range errors which are currently triggered by reading
  // the array with a concurrent write - currently Trilinos RangeError and
  // STL std::out_of_range are grouped together
  int countOutOfRangeEvents = 0;

  for (int testCycle = 0; testCycle < numTests; ++testCycle) {
    std::vector<std::thread> threads;
    ThreadTestManager::s_bAllowThreadsToRun = false;
    ThreadTestManager::s_countCompletedThreads = 0;
    ThreadTestManager::s_countWritingThreadCycles = 0;

    // 0 is pushing/popping, 1 is memory reading/writing.
    // The rest are the reader threads looking for troubles
    int finishWhenThisThreadCountCompletes = numThreads - 2;

    // I avoid using general Arrays as we are testing them,
    // makes debugging thread issues easier
    Cycle_Index_Tracker index_tracker[numThreads];

    // makes an array of length 1 - so we will cycle from size 1 to
    // size maxArraySize then back to 1
    RCP<Array<int>> array_rcp = rcp(new Array<int>(1, setValue));

    for (int i = 0; i < numThreads; ++i) {
      switch (i)
      {
        case 0:
          // create the insert/delete thread which constantly changes the
          // size of the array
          threads.push_back( std::thread(call_inserts_on_array, array_rcp,
            setValue, maxArraySize, finishWhenThisThreadCountCompletes) );
          break;
        case 1:
          // create the srambler thread which puts pressure on memory
          // with constant allocations and deletes. This allows us to detect
          // when we read bad memory (sometimes). Eventually we may fix this
          // and then the scamble events should no longer occur.
          threads.push_back( std::thread(scramble_memory, scrambleValue,
            finishWhenThisThreadCountCompletes) );
          break;
        default:
          // Default threads which just do read operations
          // We have two modes for const and non-const iterators.
          ++countTotalTestRuns;
          if (bUseConstVersion) {
            threads.push_back( std::thread(
              do_read_operations_on_array< const Array<int> >, array_rcp,
              setValue, scrambleValue, std::ref(index_tracker[i]),
              maxArraySize, arrayTestStyle));
          }
          else {
            threads.push_back( std::thread(
              do_read_operations_on_array< Array<int> >, array_rcp,
              setValue, scrambleValue, std::ref(index_tracker[i]),
              maxArraySize, arrayTestStyle));
          }
          break;
      }
    }

    // let the threads run
    ThreadTestManager::s_bAllowThreadsToRun = true;

    // join all threads
    for (unsigned int i = 0; i < threads.size(); ++i) {
      threads[i].join();
    }

    // now count all the events from the different threads
    for (unsigned int i = 0; i < threads.size(); ++i) {
      // check for a danglingReference event
      if (index_tracker[i].danglingReference != UNSET_CYCLE_INDEX ) {
        ++countDetectedDanglingReferences;
      }

      // Check for srambled memory events.
      // This is the flaw we don't currently have a solution for.
      if (index_tracker[i].scambledMemory != UNSET_CYCLE_INDEX ) {
        ++countScrambledMemoryEvents;
      }

      // Track out of range errors
      if (index_tracker[i].outOfRangeError != UNSET_CYCLE_INDEX ) {
        ++countOutOfRangeEvents;
      }

      // This version of the test is designed to hit the dangling reference
      // on the first cycle so we verify that happened - it provides a
      // mode of working which is deterministic.
      if (arrayTestStyle == ArrayTest_TriggerDanglingWithIterationFirstCycle
        && index_tracker[i].missedDanglingOnFirstCycle != UNSET_CYCLE_INDEX ) {
        ++countMissedFirstCycleDanglers;
      }
    }
    convenience_log_progress(testCycle, numTests);		// this is just output
  }

  // some log output based on the mode of the test
  switch (arrayTestStyle) {
    case ArrayTest_DoReadOperations:
    {
      // read tests are expected to hit out of range errors
      // We may also sometimes read bad memory - this currently does not
      // have a solution. It happened because the weak rcp checked it was valid,
      // thought everything was ok, then before actually reading the RCP
      // released and the scramble thread was able to change the memory.
      // If the scramble thread did not change the value, we may be reading
      // bad memory which happens to still have the right stuff.
      // So the scamble count exists to monitor that the event is taking place.
      std::cout << "Range Errors: " << countOutOfRangeEvents <<
        " Scrambles: " << countScrambledMemoryEvents << " ";
    }
    break;
    case ArrayTest_TriggerDanglingWithIterationFirstCycle:
    case ArrayTest_TriggerDanglingWithIteration:
    {
      // These two modes are expected to hit dangling references
      // and only work with the iterators.
      std::cout << "Danglers: " << countDetectedDanglingReferences << " ";
    }
    break;
  }

  // Now decide if we passed the test. For the modes with iterators,
  // meaning arrayTestStyle != ArrayTest_DoReadOperations,
  // we should have detected one dangling reference on every test run exactly.
  bool bPassed_DetectDanglers = (arrayTestStyle != ArrayTest_DoReadOperations)
    ? (countDetectedDanglingReferences == countTotalTestRuns) : true;

  // If we rang the ArrayTest_TriggerDanglingWithIterationFirstCycle mode we
  // tracked any time we did not hit the dangler on the first cycle by
  // incrementing countMissedFirstCycleDanglers. That should never happen
  // so countMissedFirstCycleDanglers should always be 0.
  bool bPassed_DetectDanglersFirstCycle = (countMissedFirstCycleDanglers == 0);

  // ArrayTest_DoReadOperations should find out of range errors so we test
  // that at least one was detected. The other two types use iterators and
  // thergore should never detect out of range errors.
  bool bPassed_CountOutOfRangeErrors =
    (arrayTestStyle == ArrayTest_DoReadOperations) ?
      (countOutOfRangeEvents != 0) : (countOutOfRangeEvents == 0);

  // however - there is an important difference for ArrayTest_DoReadOperations
  // that is the only test here which has stochastic results
  // if we are running only a few loops for the basic testing we won't
  // necessarily pick up one of those errors and then the test can fail
  // So here we should check if the test count is sufficiently high
  if(NUM_TESTS_TO_RUN < 1000 && arrayTestStyle == ArrayTest_DoReadOperations) {
    if(countOutOfRangeEvents == 0) {
      bPassed_CountOutOfRangeErrors = true; // we allow it to go through
    }
  }

  // Using the above 3 bools which were set based on the test results, we can
  // now determine if the test overall has passed. We must pass all 3.
  bool bPass = bPassed_DetectDanglersFirstCycle &&
    bPassed_DetectDanglersFirstCycle && bPassed_CountOutOfRangeErrors;

  if (!bPass) {
    std::cout << std::endl; // cosmetic - get these errors on a new line
  }

  // If we failed any of the 3 bools, we will no log some information.
  if (!bPassed_DetectDanglers) {
    std::cout << "Test FAILED because it detected only " <<
      countDetectedDanglingReferences <<
      " Dangling References but should have detected " << countTotalTestRuns
      << "." << std::endl;
  }

  if( !bPassed_DetectDanglersFirstCycle ) {
    std::cout << "Test FAILED because it missed " <<
      countMissedFirstCycleDanglers <<
      " Dangling References but should have detected " << countTotalTestRuns
      << " on the first cycle." << std::endl;
  }

  if( !bPassed_CountOutOfRangeErrors ) {
    std::cout << "Test FAILED because it detected " <<
      countOutOfRangeEvents <<
      " out of range events but should have detected: "
      << ( (arrayTestStyle == ArrayTest_DoReadOperations) ?
      "More Than 0" : "0" ) << std::endl;
  }

  return bPass;
}

// RCP Thread Safety Unit Test: mtArrayDanglingReference_NonConst_ReadValues
//
// Purpose:
//   This mode will trigger out of range errors when we might prefer to think
//   of them as dangling reference errors. This is just due to the internal
//   setup and how Array<T>::assertIndex works.
//   This may also trigger a scrambled memory event which does not have a
//   solution implemented. This was preserved to show the event is taking
//   place and future development may eliminate these events so it is
//   impossible to have a weak rcp read bad memory without detection.
//
// Description:
//   An Array is shared to multiple threads - they simultaneously try to read
//   the Array data while one thread is constantly changing the array and
//   another thread is constantly allocating and deleting memory.
//
// Solution to the Problem:
//   Same as mtArrayMultipleReads_NonConst
//
// Demonstration of Original Problem
//   Same as mtArrayMultipleReads_NonConst
TEUCHOS_UNIT_TEST( Array, mtArrayDanglingReference_NonConst_ReadValues )
{
  bool bPass = false;
  try {
    bPass = runArrayDanglingReferenceTest( false,
      ArrayTest_DoReadOperations );
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
  TEST_ASSERT( bPass )
}

// RCP Thread Safety Unit Test: mtArrayDanglingReference_Const_ReadValues
//
// Purpose:
//   Similar to mtArrayDanglingReference_NonConst_ReadValues except with
//   const iterators.
//
// Description:
//   Similar to mtArrayDanglingReference_NonConst_ReadValues
//
// Solution to the Problem:
//   Same as mtArrayMultipleReads_Const
//
// Demonstration of Original Problem
//   Same as mtArrayMultipleReads_Const
TEUCHOS_UNIT_TEST( Array, mtArrayDanglingReference_Const_ReadValues )
{
  bool bPass = false;
  try {
    bPass = runArrayDanglingReferenceTest( true,
      ArrayTest_DoReadOperations );
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
  TEST_ASSERT( bPass )
}

// RCP Thread Safety Unit Test: mtArrayDanglingReference_NonConst
//
// Purpose:
//   Demonstrates thread safe detection of dangling references which occur
//   when on thread is iterating over the array while another thread is
//   changing the size of the error (which reallocates the ArrayRCP used by
//   the iterators.
//
// Description:
//   An Array is shared to multiple threads - they simultaneously try to
//   iterate over the array while other threads change the array.
//
// Solution to the Problem:
//   Same as mtArrayMultipleReads_NonConst
//
// Demonstration of Original Problem
//   Same as mtArrayMultipleReads_NonConst
TEUCHOS_UNIT_TEST( Array, mtArrayDanglingReference_NonConst )
{
  bool bPass = false;
  try {
    bPass = runArrayDanglingReferenceTest( false,
      ArrayTest_TriggerDanglingWithIteration );
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
  TEST_ASSERT( bPass )
}

// RCP Thread Safety Unit Test: mtArrayDanglingReference_Const
//
// Purpose:
//   Same as mtArrayDanglingReference_NonConst but with const iterators.
//
// Description:
//   Similar to mtArrayDanglingReference_NonConst.
//
// Solution to the Problem:
//   Same as mtArrayMultipleReads_Const
//
// Demonstration of Original Problem
//   Same as mtArrayMultipleReads_Const
TEUCHOS_UNIT_TEST( Array, mtArrayDanglingReference_Const )
{
  bool bPass = false;
  try {
    bPass = runArrayDanglingReferenceTest( true,
      ArrayTest_TriggerDanglingWithIteration );
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
  TEST_ASSERT( bPass )
}

// RCP Thread Safety Unit Test: mtArrayDanglingReferenceFirstCycle_NonConst
//
// Purpose:
//   Same as mtArrayDanglingReference_NonConst but the test is designed
//   so the dangling reference must occur on the first cycle.
//
// Description:
//   After creating an iterator, the reading thread spin locks until it
//   determines the thread responsible for changing the array size has
//   completed at least one full cycle. Then the thread continues by using
//   the iterator which is now guaranteed to be dangling.
//
// Solution to the Problem:
//   Same as mtArrayMultipleReads_NonConst
//
// Demonstration of Original Problem
//   Same as mtArrayMultipleReads_NonConst
TEUCHOS_UNIT_TEST( Array, mtArrayDanglingReferenceFirstCycle_NonConst )
{
  bool bPass = false;
  try {
    bPass = runArrayDanglingReferenceTest( false,
      ArrayTest_TriggerDanglingWithIterationFirstCycle );
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
  TEST_ASSERT( bPass )
}

// RCP Thread Safety Unit Test: mtArrayDanglingReference_Const
//
// Purpose:
//   Similar to mtArrayDanglingReferenceFirstCycle_NonConst.
//
// Description:
//   Similar to mtArrayDanglingReferenceFirstCycle_NonConst.
//
// Solution to the Problem:
//   Same as mtArrayMultipleReads_Const
//
// Demonstration of Original Problem
//   Same as mtArrayMultipleReads_Const
TEUCHOS_UNIT_TEST( Array, mtArrayDanglingReferenceFirstCycle_Const )
{
  bool bPass = false;
  try {
    bPass = runArrayDanglingReferenceTest( true,
      ArrayTest_TriggerDanglingWithIterationFirstCycle );
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
  TEST_ASSERT( bPass )
}

#endif // TEUCHOS_DEBUG

} // end namespace



