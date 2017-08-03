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

// These unit tests are used for both a Nightly version and a Basic version

#include "General_MT_UnitTests.hpp"

#include "Teuchos_Tuple.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include <vector>
#include <thread>

namespace {

using Teuchos::Tuple;
using Teuchos::RCP;
using Teuchos::tuple;
using Teuchos::rcp;

// convenience Tuple type for testing with abritrary n = 8
#define TUPLE_SIZE 8 // arbitrary
typedef Tuple<int, TUPLE_SIZE> TupleClass;

// this class is used to test Tuple in mtTupleMultipleReads below
class TupleContainingClass
{
  public:
    TupleContainingClass()
    {
      // just set some values to match the index which we will check
      for (int n = 0; n < TUPLE_SIZE; ++n) {
        myTuple[n] = n;
      }
    }

    // this method is called by the threads and returns false if it has any
    // trouble reading the values. The constructor simply sets the elements
    // to be equal to the array index and here we verify they read that way.
    // If something happened strange with the iterators or memory we would
    // hope to trigger a problem but we did not expect this to happen and
    // everything seems ok.
    bool readTheTuple()
    {
      int sanityCheckIndex = 0;
      for (TupleClass::iterator iter = myTuple.begin();
        iter < myTuple.end(); ++iter) {
        if (sanityCheckIndex != *iter) {
          // something happened
          return false;
        }
        ++sanityCheckIndex;
      }
      return true;
    }

  private:
    TupleClass myTuple; // the Tuple to be tested

};

// utility thread method used by mtTupleMultipleReads below
static void share_tuple_to_threads(RCP<TupleContainingClass> shared_tuple,
  std::atomic<int> & countErrors) {
  while (!ThreadTestManager::s_bAllowThreadsToRun) {}
  for( int n = 0; n < 1000; ++n) {
    if (!shared_tuple->readTheTuple()) {
      ++countErrors;
    }
  }
}

// RCP Thread Safety Unit Test: mtTupleMultipleReads
//
// Purpose:
//   Sanity Check: Mirrors the Array test (which would fail without mutex
//   protection) but this oen is ok because the ArrayView it uses does not
//   have any mutable behavior - it is true const.
//
// Description:
//   Creates a Tuple, shares it to multiple threads, and reads it with
//   some verification that the reading is not doing anything odd.
//
// Solution to the Problem:
//   Sanity Check
//
// Demonstration of Problem:
//   Sanity Check
TEUCHOS_UNIT_TEST( Tuple, mtTupleMultipleReads )
{
  const int numThreads = TEUCHOS_THREAD_SAFE_UNIT_TESTS_THREADS_USED;
  const int numTests = NUM_TESTS_TO_RUN;
  // threads will increment if errors occur
  std::atomic<int> countErrors(0);
  try {
    for (int testCycle = 0; testCycle < numTests; ++testCycle) {
      std::vector<std::thread> threads;
        // note the tuple constructor makes this a type TupleClass
        // which is Tuple<int,8>
      RCP<TupleContainingClass> shared_tuple_rcp = rcp(
        new TupleContainingClass());
      // create threads all sharing access to the Tuple RCP
      for (int i = 0; i < numThreads; ++i) {
        threads.push_back( std::thread(share_tuple_to_threads,
          shared_tuple_rcp, std::ref(countErrors)) );
      }
      ThreadTestManager::s_bAllowThreadsToRun = true; // let the threads run
      for (unsigned int i = 0; i < threads.size(); ++i) {
        threads[i].join();
      }
      convenience_log_progress(testCycle, numTests); // this is just output
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
  TEST_EQUALITY(countErrors, 0 );
}

} // end namespace
