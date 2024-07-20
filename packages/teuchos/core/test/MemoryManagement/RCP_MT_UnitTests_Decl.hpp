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

#include "Teuchos_RCP.hpp"
#include "Teuchos_RCPNode.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include <vector>
#include <thread>

namespace {

using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;

// Thread Utility Method
// See unit test mtRefCount below for details.
// This method is called by several threads so each can copy the same RCP
static void make_large_number_of_copies(RCP<int> ptr) {
  std::vector<RCP<int> > ptrs(10000, ptr);
}

// RCP Thread Safety Unit Test: mtRefCount
//
// Purpose:
//   Test that RCP counters are thread safe.
//
// Description:
//   This was the first unit test created for developing
//   RCP thread safe code. The purpose of the test was to demonstrate that the
//   int RCP counters for weak and strong were not thread safe. This was not
//   unexpected because they were not atomics.
//   In this test, several threads are all passed an RCP<int> and simply
//   create many copies of the same RCP. They all share the same root node
//   and the total strong count should be consistent with the number of copies.
//   As the threads complete, they release all of the RCP objects and the
//   strong counters should deincrement accordingly. At the end of the test
//   the total strong count should be 1 because there is only the original
//   RCP<int> remaining.
//   If the counters were not atomic, simultaneous ++i and --i action on the int
// counters will result in a a failed total count.
//
// Solution to the Problem:
//   This issue was resolved by changing the counters
//   to be std::atomic<int>. This is expected to have a performance impact.
//
// Demonstration of Problem:
//   To reproduce the original problem, add the following define at the
//   top of this file:
//     #define DISABLE_ATOMIC_COUNTERS
//   That wil change the counters back to  int instead of atomic<int>
//   and the test will then fail.
TEUCHOS_UNIT_TEST( RCP, mtRefCount )
{
  const int numThreads = TEUCHOS_THREAD_SAFE_UNIT_TESTS_THREADS_USED;
  RCP<int> ptr(new int); // RCP that all the threads will copy
  std::vector<std::thread> threads;
  for (int i = 0; i < numThreads; ++i) {
    threads.push_back(std::thread(make_large_number_of_copies, ptr));
  }
  for (unsigned int i = 0; i < threads.size(); ++i) {
    threads[i].join();
  }
  // we still have one strong RCP so the total_count should be 1
  TEST_EQUALITY_CONST(ptr.total_count(), 1);
}

// Thread Utility Method
// See unit test mtCreateIndependentRCP below for description.
// Each thread makes new RCP<int> objects and deletes them to stress test
// the RCPNodeTracer mechanism
static void create_independent_rcp_objects() {
  // spin lock the threads so we can trigger them all together
  while (!ThreadTestManager::s_bAllowThreadsToRun) {}
  for(int n = 0; n < NUM_TESTS_TO_RUN; ++n ) {
    // this allocates a new rcp ptr independent of all other rcp ptrs,
    // and then dumps it, over and over
    RCP<int> ptr( new int );
  }
}

// RCP Thread Safety Unit Test: mtCreateIndependentRCP
//
// Purpose:
//   Test that the debug RCPNodeTracer is thread safe.
//
// Description:
//   Multiple threads all create their own individual RCP<int> objects and
//   release them many times. In Debug mode, the RCPNodeTracer would fail.
//   Note that in release this test runs but does not test anything important.
//
// Solution to the Problem:
//   Used a static mutex in Teuchos_RCPNode.cpp
//
// Demonstration of Problem:
//   To reproduce the original problem, add the following define to the top of
//   the file: Teuchos_RCPNode.cpp
//     #define USE_MUTEX_TO_PROTECT_NODE_TRACING
//   That wil remove the mutex which was added to protect RCPNodeTracker and
//   restore the original errors. Note that will lead to undefined behaviors.
TEUCHOS_UNIT_TEST( RCP, mtCreateIndependentRCP )
{
  const int numThreads = TEUCHOS_THREAD_SAFE_UNIT_TESTS_THREADS_USED;
  // track node count when we start because other RCP obejcts already exist
  int initialNodeCount = Teuchos::RCPNodeTracer::numActiveRCPNodes();
  try {
    std::vector<std::thread> threads;
    // threads will start in a spin-lock state
    ThreadTestManager::s_bAllowThreadsToRun = false;
    for (int i = 0; i < numThreads; ++i) {
      threads.push_back(std::thread(create_independent_rcp_objects));
    }
    // releasa all the threads
    ThreadTestManager::s_bAllowThreadsToRun = true;
    for (unsigned int i = 0; i < threads.size(); ++i) {
      threads[i].join();
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
  // we should be back to where we were when we started the test
  // this is not going to be 0 because other objects existed before the test
  TEST_EQUALITY_CONST(initialNodeCount,
    Teuchos::RCPNodeTracer::numActiveRCPNodes());
}

// Thread Utility Method
// See unit test mtTestGetExistingRCPNodeGivenLookupKey below for details.
static void create_independent_rcp_without_ownership() {
  for(int n = 0; n < 10000; ++n ) {
    int * intPtr = new int;
    // this allocates a new rcp but without memory ownership
    RCP<int> ptr( intPtr, false );
    // since the RCP doesn't have memory ownership, we delete it manually
    delete intPtr;
  }
}

// RCP Thread Safety Unit Test: mtTestGetExistingRCPNodeGivenLookupKey
//
// Purpose:
//   Test that getExistingRCPNodeGivenLookupKey() is thread safe.
//
// Description:
//   Multiple threads all create their own individual RCP<int> objects
//   without ownership and then handle deleting the memory.
//
// Solution to the Problem:
//   Added mutex protection to getExistingRCPNodeGivenLookupKey()
//
// Demonstration of Problem:
//   Comment out the lock_guard in getExistingRCPNodeGivenLookupKey() in
//   Teuchos_RCPNode.cpp which will lead to undefined behaviors during the test.
TEUCHOS_UNIT_TEST( RCP, mtTestGetExistingRCPNodeGivenLookupKey )
{
  const int numThreads = TEUCHOS_THREAD_SAFE_UNIT_TESTS_THREADS_USED;
  try {
    std::vector<std::thread> threads;
    for (int i = 0; i < numThreads; ++i) {
      threads.push_back(std::thread(create_independent_rcp_without_ownership));
    }
    for (unsigned int i = 0; i < threads.size(); ++i) {
      threads[i].join();
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
}

// Thread Utility Method
// This is used by several of the unit tests below.
// It simply provides an RCP to a thread.
// When the thread is destroyed, it releases the passes RCP object.
// See the test descriptions for more details.
template<typename SOURCE_RCP_TYPE>
static void thread_gets_a_copy_of_rcp(SOURCE_RCP_TYPE ptr) {
  // spin lock the threads so we can trigger them all at once
  // note we don't actually do anything - the thread was passed a copy which is
  // all we need for this test - it will be deleted when the thread ends.
  while(!ThreadTestManager::s_bAllowThreadsToRun) {}
}

// RCP Thread Safety Unit Test: mtRCPLastReleaseByAThread
//
// Purpose:
//   Demonstrate that the reading of the deincr_count() in RCP is thread safe.
//
// Description:
//   Demonstrate the delete path was not thread safe - specifically had to
//   modify the format of deincr_count(). Calling --count and then checking
//   count was not thread safe so we must check if(--count ==0) to have a
//   thread safe atomic read of the deincrement event.
//
// Solution to the Problem:
//   Changed the logic of the deincrement() counter path to be atomically safe.
//
// Demonstration of Problem:
//   Add the following define to the top of this file:
//     #define BREAK_THREAD_SAFETY_OF_DEINCR_COUNT
//   That will create an error in the RCP code which mimics the way the code
//   originally ran before the development of the thread safe system.
//   Note that will likely lead to undefined behaviors in removeRCPNode().
TEUCHOS_UNIT_TEST( RCP, mtRCPLastReleaseByAThread )
{
  const int numThreads = TEUCHOS_THREAD_SAFE_UNIT_TESTS_THREADS_USED;
  const int numCycles = NUM_TESTS_TO_RUN;
  try {
    for(int cycleIndex = 0; cycleIndex < numCycles; ++cycleIndex) {
      // initialize
      CatchMemoryLeak::s_countAllocated = 0;
      // only 1 new allocation happens in this test
      RCP<CatchMemoryLeak> ptr(new CatchMemoryLeak);
      // prepare to spin lock the threads
      ThreadTestManager::s_bAllowThreadsToRun = false;
      std::vector<std::thread> threads;
      for (int threadIndex = 0; threadIndex < numThreads; ++threadIndex) {
        threads.push_back(std::thread(
          thread_gets_a_copy_of_rcp<RCP<CatchMemoryLeak>>, ptr));
      }
      // at this point threads are spin locked and holding copies
      // release the ptr in the main thread
      ptr = null;
      // now we release all the threads
      ThreadTestManager::s_bAllowThreadsToRun = true;
      for (unsigned int i = 0; i < threads.size(); ++i) {
        // when join completes rcp should be completely deleted
        threads[i].join();
      }
      convenience_log_progress(cycleIndex, numCycles);	// this is just output
      if (CatchMemoryLeak::s_countAllocated != 0) {
        break; // will catch in error below
      }
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
  // test for valid RCP deletion
  TEST_EQUALITY_CONST(CatchMemoryLeak::s_countAllocated, 0);
}

// This dealloc method is passed to RCP objects to be called when they delete.
// To track the behavior an atomic counter records all the events so that
// thread issues can be detected.
// This method is used by unit test mtRCPLastReleaseByAThreadWithDealloc below.
void deallocCatchMemoryLeak(CatchMemoryLeak* ptr)
{
  // increment the static global counter used by the unit tests
  ++CatchMemoryLeak::s_countDeallocs;
  // then implement the delete of the data
  delete ptr;
}

// RCP Thread Safety Unit Test: mtRCPLastReleaseByAThreadWithDealloc
//
// Purpose:
//   Sanity Check: Similar to mtRCPLastReleaseByAThread
//
// Description:
//   This unit test is identical to mtRCPLastReleaseByAThread except that the
//   RCP now has a custom dealloc policy. Once the delete path was made thread
//   safe, this is also expected to be thread safe, so no errors were expected
//   and no additional probelms were observed.
//
// Solution to the Problem:
//   Sanity Check
//
// Demonstration of Problem:
//   Sanity Check
TEUCHOS_UNIT_TEST( RCP, mtRCPLastReleaseByAThreadWithDealloc )
{
  const int numThreads = TEUCHOS_THREAD_SAFE_UNIT_TESTS_THREADS_USED;
  const int numCycles = NUM_TESTS_TO_RUN;
  try {
    for(int cycleIndex = 0; cycleIndex < numCycles; ++cycleIndex) {
      CatchMemoryLeak::s_countDeallocs = 0; // set it to 0
      RCP<CatchMemoryLeak> ptr = rcpWithDealloc(new CatchMemoryLeak,
        Teuchos::deallocFunctorDelete<CatchMemoryLeak>(deallocCatchMemoryLeak));
      // prepare to spin lock the threads
      ThreadTestManager::s_bAllowThreadsToRun = false;
      std::vector<std::thread> threads;
      for (int threadIndex = 0; threadIndex < numThreads; ++threadIndex) {
        threads.push_back(std::thread(
          thread_gets_a_copy_of_rcp<RCP<CatchMemoryLeak>>, ptr));
      }
      // at this point threads are spin locked and holding copies
      // Release the ptr in the main thread
      ptr = null;
      // now we release all the threads
      ThreadTestManager::s_bAllowThreadsToRun = true;
      for (unsigned int i = 0; i < threads.size(); ++i) {
        // when join completes rcp should be completely deleted
        threads[i].join();
      }
      convenience_log_progress(cycleIndex, numCycles);// this is just output
      if (CatchMemoryLeak::s_countDeallocs != 1) {
        break; // will catch in error below
      }
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
  // we should have ended with exactly one dealloc call
  TEST_EQUALITY_CONST(CatchMemoryLeak::s_countDeallocs, 1);
}

// This method is used by mtRCPLastReleaseByAThreadWithDeallocHandle below.
void deallocHandleCatchMemoryLeak(CatchMemoryLeak** handle)
{
  ++CatchMemoryLeak::s_countDeallocs;
  CatchMemoryLeak *ptr = *handle;
  delete ptr;
  *handle = 0;
}


// RCP Thread Safety Unit Test: mtRCPLastReleaseByAThreadWithDeallocHandle
//
// Purpose:
//   Sanity Check: Similar to mtRCPLastReleaseByAThread
//
// Description:
//   This unit test is identical to mtRCPLastReleaseByAThread except that the
//   RCP now has a custom dealloc handle policy. Similar to the test
//   mtRCPLastReleaseByAThreadWithDealloc above, this was not found to have
//   any additional problems.
//
// Solution to the Problem:
//   Sanity Check
//
// Demonstration of Problem:
//   Sanity Check
TEUCHOS_UNIT_TEST( RCP, mtRCPLastReleaseByAThreadWithDeallocHandle )
{
  const int numThreads = TEUCHOS_THREAD_SAFE_UNIT_TESTS_THREADS_USED;
  const int numCycles = NUM_TESTS_TO_RUN;
  try {
    for(int cycleIndex = 0; cycleIndex < numCycles; ++cycleIndex) {
      CatchMemoryLeak::s_countDeallocs = 0; // set it to 0
      RCP<CatchMemoryLeak> ptr = rcpWithDealloc(new CatchMemoryLeak,
       Teuchos::deallocFunctorHandleDelete<CatchMemoryLeak>(
         deallocHandleCatchMemoryLeak));
      // prepare the threads to be spin locked
      ThreadTestManager::s_bAllowThreadsToRun = false;
      std::vector<std::thread> threads;
      for (unsigned int threadIndex = 0; threadIndex <
        numThreads; ++threadIndex) {
        threads.push_back(std::thread(
          thread_gets_a_copy_of_rcp<RCP<CatchMemoryLeak>>, ptr));
      }
      // at this point threads are spin locked and holding copies
      // Release the ptr in the main thread
      ptr = null;
      // now we release all the threads
      ThreadTestManager::s_bAllowThreadsToRun = true;
      for (unsigned int i = 0; i < threads.size(); ++i) {
        // when join completes rcp should be completely deleted
        threads[i].join();
      }
      convenience_log_progress(cycleIndex, numCycles); // this is just output
      if (CatchMemoryLeak::s_countDeallocs != 1) {
        break; // will catch in error below
      }
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
  TEST_EQUALITY_CONST(CatchMemoryLeak::s_countDeallocs, 1);	// should be 1 only
}


// See mtRCPThreadCallsRelease below for details
// This utitily method allows one thread to call release on the RCP
static void call_release_on_rcp_if_flag_is_set(RCP<CatchMemoryLeak> ptr,
  int numCopies, bool bCallsRelease) {
  // spin lock the threads so we can trigger them all at once
  while(!ThreadTestManager::s_bAllowThreadsToRun) {}
  if(bCallsRelease) {
    ptr.release(); // should make a memory leak!
  }
}

// RCP Thread Safety Unit Test: mtRCPThreadCallsRelease
//
// Purpose:
//   Sanity Check: To validate the release mechanism. There was no explicit
//   problem expected and none was observed.
//
// Description:
//   In this test only one thread calls release which means the data should
//   not be deleted. It is expected to be deleted manually. The test verifies
//
// Solution to the Problem:
//   Sanity Check
//
// Demonstration of Problem:
//   Sanity Check
TEUCHOS_UNIT_TEST( RCP, mtRCPThreadCallsRelease )
{
  const int numThreads = TEUCHOS_THREAD_SAFE_UNIT_TESTS_THREADS_USED;
  const int numCycles = NUM_TESTS_TO_RUN;
  bool bFailure = false;
  try {
    for(int cycleIndex = 0; cycleIndex < numCycles; ++cycleIndex) {
      CatchMemoryLeak::s_countAllocated = 0; // initialize
      CatchMemoryLeak * pMemoryToLeak = new CatchMemoryLeak;
      RCP<CatchMemoryLeak> ptr(pMemoryToLeak);
      // prepare to spin lock the threads
      ThreadTestManager::s_bAllowThreadsToRun = false;
      std::vector<std::thread> threads;
      for (int threadIndex = 0; threadIndex < numThreads; ++threadIndex) {
        bool bCallRelease = (threadIndex==0); // only the first calls release
        threads.push_back(std::thread(call_release_on_rcp_if_flag_is_set,
          ptr, 1, bCallRelease));
      }
      ptr = null;
      ThreadTestManager::s_bAllowThreadsToRun = true;
      for (unsigned int i = 0; i < threads.size(); ++i) {
        threads[i].join();
      }
      convenience_log_progress(cycleIndex, numCycles); // this is just output
      if (CatchMemoryLeak::s_countAllocated != 1) {
        break; // will catch in error below
      }
      else {
        // we can clean up our memory leak now by hand
        // if the test is passing we want a clean finish
        delete pMemoryToLeak;
      }
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
  // organized like this because I want the above loops to properly clean up
  // memory if the test is passing - which changes the counters and makes this
  // interpretation complicated
  TEST_EQUALITY_CONST(bFailure, false);
}

// Used by unit test mtRCPExtraData below to test ExtraData feature of RCP
template<typename T>
class ExtraDataTest {
public:
  static RCP<ExtraDataTest<T> > create(T *ptr)
    { return rcp(new ExtraDataTest(ptr)); }
  ~ExtraDataTest() { delete [] ptr_; } // proper delete
private:
  T *ptr_;
  ExtraDataTest(T *ptr) : ptr_(ptr) {}
  // Not defined!
  ExtraDataTest();
  ExtraDataTest(const ExtraDataTest&);
  ExtraDataTest& operator=(const ExtraDataTest&);
};

// RCP Thread Safety Unit Test: mtRCPExtraData
//
// Purpose:
//   Sanity Check: Use the extra data feature and make sure things are ok.
//   There was no explicit problem expected and none was observed.
//
// Description:
//   In this test the ExtraData RCP feature is used to properly corred the ptr
//   to delete using delete [].
//
// Solution to the Problem:
//   Sanity Check
//
// Demonstration of Problem:
//   Sanity Check
TEUCHOS_UNIT_TEST( RCP, mtRCPExtraData )
{
  const int numThreads = TEUCHOS_THREAD_SAFE_UNIT_TESTS_THREADS_USED;
  const int numCycles = NUM_TESTS_TO_RUN;
  try {
    for(int cycleIndex = 0; cycleIndex < numCycles; ++cycleIndex) {
      CatchMemoryLeak::s_countAllocated = 0; // initialize
      // standard delete will be wrong - should call delete[]
      RCP<CatchMemoryLeak> ptr(new CatchMemoryLeak[1]);
      ptr.release();		// extra data will handle the case
      Teuchos::set_extra_data( ExtraDataTest<CatchMemoryLeak>::create(
        ptr.getRawPtr()), "dealloc", Teuchos::inOutArg(ptr));
        // prepare to spin lock the threads
      ThreadTestManager::s_bAllowThreadsToRun = false;
      std::vector<std::thread> threads;
      for (unsigned int threadIndex = 0; threadIndex < numThreads;
        ++threadIndex) {
        threads.push_back(std::thread(
          thread_gets_a_copy_of_rcp<RCP<CatchMemoryLeak>>, ptr));
      }
      // At this point threads are spin locked and holding copies
      // Release the ptr in the main thread
      ptr = null;
      // now we release all the threads
      ThreadTestManager::s_bAllowThreadsToRun = true;
      for (unsigned int i = 0; i < threads.size(); ++i) {
        // when join completes rcp should be completely deleted
        threads[i].join();
      }
      convenience_log_progress(cycleIndex, numCycles); // this is just output
      if (CatchMemoryLeak::s_countAllocated != 0) {
        break; // will catch in error below
      }
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
  // test for valid RCP deletion
  TEST_EQUALITY_CONST(CatchMemoryLeak::s_countAllocated, 0);
}

// RCP Thread Safety Unit Test: mtRCPWeakStrongDeleteRace
//
// Purpose:
//   To demonstrate weak and strong RCP objects can delete simultaneously
//   and be handled in a thread safe way.
//
// Description:
//   4 threads are creating alternating between having a strong and weak RCP.
//   They all share the same root RCP object.
//   The main thread will release the RCP and then release all the threads.
//   When the threads die the RCP objects will release and strong/weak
//   counters will process in a race condition enivronment.
//
// Solution to the Problem:
//   counters were refactored to work like boost, which uses a nice trick where
//   the strong count is as expected, but the weak count is equal to weak + 1.
//   This allows all the race conditions to be resolved in an elegant way.
//
// Demonstration of Problem:
//   I had a define for this but it seemed really messy to leave that in core
//   So here is a way to simulate the original problem:
//   In Teuchos_RCPNode.hpp change the unbind() method:
//   Change that one line unbindOneStrong() to be:
//            node_->deincr_count(RCP_WEAK);
//            unbindOneStrong();
//            node_->incr_count(RCP_WEAK);
//   That will give a behavior similar to the original code.
TEUCHOS_UNIT_TEST( RCP, mtRCPWeakStrongDeleteRace )
{
  const int numThreads = TEUCHOS_THREAD_SAFE_UNIT_TESTS_THREADS_USED;
  const int numCycles = NUM_TESTS_TO_RUN;
  try {
    for(int cycleIndex = 0; cycleIndex < numCycles; ++cycleIndex) {
      CatchMemoryLeak::s_countAllocated = 0; // initialize
      // only 1 new allocation happens in this test
      RCP<CatchMemoryLeak> ptr(new CatchMemoryLeak);
      // prepare to spin lock the threads
      ThreadTestManager::s_bAllowThreadsToRun = false;
      std::vector<std::thread> threads;
      bool bToggleStrong = true;
      for (int threadIndex = 0; threadIndex < numThreads; ++threadIndex) {
        if (bToggleStrong) {
          threads.push_back(std::thread(
            thread_gets_a_copy_of_rcp<RCP<CatchMemoryLeak>>,
              ptr.create_strong()));
        }
        else {
          threads.push_back(std::thread(
            thread_gets_a_copy_of_rcp<RCP<CatchMemoryLeak>>,
              ptr.create_weak()));
        }
        bToggleStrong = !bToggleStrong;
      }
      // at this point threads are spin locked and holding copies
      // Release the ptr in the main thread
      ptr = null;
      // now we release all the threads
      ThreadTestManager::s_bAllowThreadsToRun = true;
      for (unsigned int i = 0; i < threads.size(); ++i) {
        // when join completes rcp should be completely deleted
        threads[i].join();
      }
      convenience_log_progress(cycleIndex, numCycles);	// this is just output
      if (CatchMemoryLeak::s_countAllocated != 0) {
        break; // will catch in error below
      }
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
  // test for valid RCP deletion
  TEST_EQUALITY_CONST(CatchMemoryLeak::s_countAllocated, 0);
}

// This utlity method supports the unit test mtRCPWeakStrongDeleteRace below.
static std::atomic<int> s_count_successful_conversions(0);
static std::atomic<int> s_count_failed_conversions(0);
template<class SOURCE_RCP_TYPE>
static void attempt_make_a_strong_ptr(SOURCE_RCP_TYPE ptr) {
  // spin lock the threads so we can trigger them all at once
  while(!ThreadTestManager::s_bAllowThreadsToRun) {}
  // ptr can be weak or strong - the weak ptrs may fail
  RCP<CatchMemoryLeak> possibleStrongPtr = ptr.create_strong_thread_safe();
  if (possibleStrongPtr.is_null()) {
    ++s_count_failed_conversions;
  }
  else {
    ++s_count_successful_conversions;
  }
}

// RCP Thread Safety Unit Test: mtRCPWeakStrongDeleteRace
//
// Purpose:
//   To demonstrate thread safe conversion of weak to strong ptr.
//   NOTE: This is not currently part of the primary goals - it is expected
//   to have application in future work.
//
// Description:
//   The threads alternate to have weak or strong ptrs to the same RCP.
//   The main thread releases it's RCP.
//   Now all the thread attempt to make a strong ptr and release their threads.
//   This creates a mixed race condition where the attempt to create a strong
//   ptr may fail if all other strong RCP's happen to be gone. This should
//   handle in a thread safe way and the failure should give a clean null result.
//
// Solution to the Problem:
//   attemptIncrementStrongCountFromNonZeroValue() was created to implement
//   this in similar fashion to the boost methods.
//
// Demonstration of Problem:
//   To see the test fail, we can modify the Teuchos_RCPNode.hpp method
//   attemptIncrementStrongCountFromNonZeroValue() to behave as if USING_ATOMICS
//   was not defined. Then the test will detect the failures.
TEUCHOS_UNIT_TEST( RCP, mtRCPMixedWeakAndStrongConvertToStrong )
{
  const int numThreads = TEUCHOS_THREAD_SAFE_UNIT_TESTS_THREADS_USED;
  int numCycles = NUM_TESTS_TO_RUN;
  s_count_successful_conversions = 0;
  s_count_failed_conversions = 0;
  try {
    for(int cycleIndex = 0; cycleIndex < numCycles; ++cycleIndex) {
      CatchMemoryLeak::s_countAllocated = 0; // initialize
      // only 1 new allocation happens in this test
      RCP<CatchMemoryLeak> ptr(new CatchMemoryLeak);
      // prepare to spin lock the threads
      ThreadTestManager::s_bAllowThreadsToRun = false;
      std::vector<std::thread> threads;
      bool bCycleStrong = true;
      for (int threadIndex = 0; threadIndex < numThreads; ++threadIndex) {
        if (bCycleStrong) {
          threads.push_back(std::thread(
            attempt_make_a_strong_ptr<RCP<CatchMemoryLeak>>,
              ptr.create_strong()));
        }
        else {
          threads.push_back(std::thread(
            attempt_make_a_strong_ptr<RCP<CatchMemoryLeak>>,
              ptr.create_weak()));
        }
        bCycleStrong = !bCycleStrong;
      }
      // at this point threads are spin locked and holding copies
      // Release the ptr in the main thread
      ptr = null;
      // now we release all the threads
      ThreadTestManager::s_bAllowThreadsToRun = true;
      for (unsigned int i = 0; i < threads.size(); ++i) {
        // when join completes rcp should be completely deleted
        threads[i].join();
      }
      if (CatchMemoryLeak::s_countAllocated != 0) {
        break;
      }
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
  std::cout << "Weak converted with null " << s_count_failed_conversions <<
    " times and success " << s_count_successful_conversions
    << " times. We want to see a mix of each. ";

  // for nightly this should definitely hit both fails and success
  // note that here 'failed' does not mean bad - it means the test properly
  // determined that the weak to strong conversion was not possible.
  // for basic I am disabling this check because we just run a few times and
  // it's possible we won't get enough checks to be sure to hit each type
  if(NUM_TESTS_TO_RUN >= 100) {
    // this has to be a mixed result or the test is not doing anything useful
    TEST_INEQUALITY_CONST(s_count_failed_conversions, 0);
    // this has to be a mixed result or the test is not doing anything useful
    TEST_INEQUALITY_CONST(s_count_successful_conversions, 0);
  }

  TEST_EQUALITY(CatchMemoryLeak::s_countAllocated, 0);     // should be 0
}

} // namespace
