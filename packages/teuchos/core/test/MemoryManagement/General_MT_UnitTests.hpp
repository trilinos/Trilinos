// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_GENERAL_MT_UNITTESTS_HPP
#define TEUCHOS_GENERAL_MT_UNITTESTS_HPP

#include "TeuchosCore_ConfigDefs.hpp"
#include "Teuchos_ConfigDefs.hpp"

// NUM_TOTAL_CORES_USED = 4 set in CMake
#define TEUCHOS_THREAD_SAFE_UNIT_TESTS_THREADS_USED 4

#include <atomic>
#include <iostream>

namespace {

// this is a convenience function for outputting percent complete information
// for long tests designed to find race conditions
static void convenience_log_progress(int cycle, int totalCycles) {
  if (cycle==0) {
    std::cout << "Percent complete: "; // begin the log line
  }
  // log every 10% percent complete - using mod % to output at regular intervals
  int mod = (totalCycles/10);
  // sometimes quick testing so make sure mod is not 0
  if((cycle % (mod == 0 ? 1 : mod) == 0) || (cycle == totalCycles-1)) {
    std::cout <<
      (int)( 100.0f * (float) cycle / (float) (totalCycles-1) ) << "% ";
    // not necessary on some setups but for Xcode this would not work without
    // the flush command - it waits for an endl
    std::flush( std::cout );
  }
}

// This class is a utility class which tracks constructor/destructor calls
// (for this test) or counts times a dealloc or deallocHandle was implemented
// (for later tests)
class CatchMemoryLeak
{
public:
  CatchMemoryLeak() { ++s_countAllocated; }
  ~CatchMemoryLeak() { --s_countAllocated; }
  static std::atomic<int> s_countAllocated;
  static std::atomic<int> s_countDeallocs;
};
// counts constructor calls (+1) and destructor calls (-1) which may include
// double delete events
std::atomic<int> CatchMemoryLeak::s_countAllocated(0);
// counts dealloc or dellocHandle calls - used for test 4 and test 5
std::atomic<int> CatchMemoryLeak::s_countDeallocs(0);

// manages a bool for spin locking unit test threads - the idea is to hold all
// threads until ready and then release simultaneously to maximize race
// condition probability.
class ThreadTestManager
{
public:
  // used to spin lock the sub threads until all released simultaneously
  static std::atomic<bool> s_bAllowThreadsToRun;
  // this lets the sub threads know the event happened so they can quit or give
  // errors if they don't detect something after this occurs
  static std::atomic<bool> s_bMainThreadSetToNull;
  // utility to tell special threads when main threads have finished work
  static std::atomic<int> s_countCompletedThreads;
  // utility to count how many times a working thread has run - using this to
  // debug other threads if they are supposed to trigger an event once this
  // thread completes a cycle
  static std::atomic<int> s_countWritingThreadCycles;
};
std::atomic<bool> ThreadTestManager::s_bAllowThreadsToRun(false);
std::atomic<bool> ThreadTestManager::s_bMainThreadSetToNull(false);
std::atomic<int> ThreadTestManager::s_countCompletedThreads(0);
std::atomic<int> ThreadTestManager::s_countWritingThreadCycles(0);

// utility define used below: -1 means it was never set
#define UNSET_CYCLE_INDEX -1

// this convenience struct is used to track the various errors which can
// happen with the RCP classes, which uses weak references. The unit test
// description below describes in more detail what this test is doing.
struct Cycle_Index_Tracker
{
  Cycle_Index_Tracker()
  {
    missedDanglingOnFirstCycle = UNSET_CYCLE_INDEX;
    danglingReference = UNSET_CYCLE_INDEX;
    scambledMemory = UNSET_CYCLE_INDEX;
    outOfRangeError = UNSET_CYCLE_INDEX;
    unknownError = UNSET_CYCLE_INDEX;
  }
  // for special test with detection designed to hit on the first cycle,
  // this tracks if it actually missed that cycle
  // this test has less random behavior so was useful for validation
  int missedDanglingOnFirstCycle;
  // tracks when a dangling reference was hit
  int danglingReference;
  // tracks when scrambled memory was detected
  int scambledMemory;
  // tracks when an out of range error was found
  int outOfRangeError;
  // tracks unknown exception - not expected
  int unknownError;
  // feedback to indicate how many times the thread has actually done a loop
  std::atomic<int> trackCycle;
};

} // end namespace

#endif // end #ifdef TEUCHOS_GENERAL_MT_UNITTESTS_HPP
