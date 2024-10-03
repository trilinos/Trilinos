// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_TimeMonitor.hpp"

namespace { // (anonymous)

using Teuchos::Array;
using Teuchos::outArg;
using Teuchos::RCP;
using Teuchos::REDUCE_MAX;
using Teuchos::reduceAll;
using std::cerr;
using std::cout;
using std::endl;

// This test works (and exercises the interesting case) in serial mode
// or for 1 MPI process, but it was originally written for 2 MPI
// processes.
TEUCHOS_UNIT_TEST( Map, Bug5822_StartWithZeroThenSkipTo3Billion )
{
  int locallyThrew = 0;
  int globallyThrew = 0;

  auto comm = Tpetra::getDefaultComm ();
  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();
  TEUCHOS_TEST_FOR_EXCEPTION(numProcs != 2, std::logic_error,
    "This test only makes sense to run with 2 MPI processes.");

  // Pick the global index type to have at least 64 bits.
#if ! defined (HAVE_TPETRA_INT_LONG_LONG) && ! defined (HAVE_TPETRA_INT_LONG)
  using GO = Tpetra::Map<>::global_ordinal_type; // just to make it compile
  out << "This test only makes sense if GO = long or long long is enabled.  "
    "That is because the test is supposed to exercise global indices greater "
    "than the maximum that int can represent (about 2 billion)." << endl;
  return;
#else
#  if defined (HAVE_TPETRA_INT_LONG_LONG)
  using GO = long long; // C++11 requires that sizeof(long long) >= 8
#  elif defined (HAVE_TPETRA_INT_LONG)
  using GO = long; // long is 32 bits on some platforms, including Windows
  if (sizeof (long) <= 4) {
    out << "sizeof (long) = " << sizeof (long) << " <= 4.  "
      "This test only makes sense if sizeof (long) >= 8.  "
      "That is because the test is supposed to exercise global indices "
      "greater than the maximum that int can represent (about 2 billion)."
        << endl;
    return;
  }
#  endif
#endif
  using LO = Tpetra::Map<>::local_ordinal_type;
  using map_type = Tpetra::Map<LO, GO>;

  // Proc 0 gets [0, 3B, 3B+2, 3B+4, 3B+6, 3B+8] (6 GIDs).
  // Proc 1 gets [3B+12, 3B+14, 3B+16, 3B+18, 3B+20] (5 GIDs).
  //
  // The GIDs are not contiguous in order to prevent Map from
  // detecting contiguity and bypassing the noncontiguous Map
  // case.
  const size_t localNumElts = (myRank == 0) ? 6 : 5;
  const Tpetra::global_size_t globalNumElts = 1 + 5*comm->getSize ();
  const GO globalFirstGid = 1L;
  const GO threeBillion = static_cast<GO> (3000000000L);

  Array<GO> myGids (localNumElts);
  // Make a copy, just to make sure that Map's constructor didn't
  // sneakily change the input ArrayView.
  Array<GO> myGidsExpected (localNumElts);
  if (myRank == 0) {
    myGids[0] = globalFirstGid;
    myGidsExpected[0] = myGids[0];
    myGids[1] = threeBillion;
    myGidsExpected[1] = myGids[1];
    for (size_t k = 2; k < localNumElts; ++k) {
      myGids[k] = myGids[k-1] + 2;
      myGidsExpected[k] = myGids[k];
    }
  } else {
    myGids[0] = threeBillion + Teuchos::as<GO> ((localNumElts+1) * 2);
    myGidsExpected[0] = myGids[0];
    for (size_t k = 1; k < localNumElts; ++k) {
      myGids[k] = myGids[k-1] + 2;
      myGidsExpected[k] = myGids[k];
    }
  }

  // if (myRank == 0) {
  //   std::cout << "type '0' and hit enter" << std::endl;
  //   int zero;
  //   std::cin >> zero;
  // }
  // comm->barrier();
  // Tpetra::Map requires that the index base is less than the minimum GID.
  const GO indexBase = 0L;
  out << "Constructing the Map" << endl;
  RCP<const map_type> map;
  try {
    TEUCHOS_FUNC_TIME_MONITOR("Construct Map");
    map = Teuchos::rcp (new map_type (globalNumElts, myGids (), indexBase, comm));
  }
  catch (std::exception& e) {
    locallyThrew = 1;
    std::ostringstream os;
    os << "Proc " << myRank << ": At Map creation, locally threw exception: "
       << e.what () << endl;
    cerr << os.str ();
  }
  globallyThrew = 0;
  reduceAll<int, int> (*comm, REDUCE_MAX, locallyThrew, outArg (globallyThrew));
  TEST_EQUALITY_CONST( globallyThrew, 0 );

  cerr << myRank << ": Querying the Map for local elements" << endl;
  {
    TEUCHOS_FUNC_TIME_MONITOR("Querying the Map for local elements");
    Teuchos::ArrayView<const GO> myGidsFound = map->getLocalElementList ();
    TEST_COMPARE_ARRAYS( myGidsExpected (), myGidsFound () );
  }

  cerr << myRank << ": Querying the Map for remote elements" << endl;
  // Proc 0 gets [1, 3B, 3B+2, 3B+4, 3B+6, 3B+8] (6 GIDs).
  // Proc 1 gets [3B+12, 3B+14, 3B+16, 3B+18, 3B+20] (5 GIDs).
  {
    TEUCHOS_FUNC_TIME_MONITOR("Querying the Map for remote elements");

    const int numRemoteGids = (myRank == 0) ? 5 : 6;

    Array<GO> remoteGids (numRemoteGids);
    Array<int> remotePids (numRemoteGids, -1);
    Array<LO> remoteLids (numRemoteGids, Teuchos::OrdinalTraits<LO>::invalid ());
    if (myRank == 0) {
      try {
        Array<int> expectedRemotePids (numRemoteGids);
        std::fill (expectedRemotePids.begin (), expectedRemotePids.end (), 1);
        Array<int> expectedRemoteLids (numRemoteGids);
        expectedRemoteLids[0] = 0;
        expectedRemoteLids[1] = 1;
        expectedRemoteLids[2] = 2;
        expectedRemoteLids[3] = 3;
        expectedRemoteLids[4] = 4;
        remoteGids[0] = threeBillion + 12;
        remoteGids[1] = threeBillion + 14;
        remoteGids[2] = threeBillion + 16;
        remoteGids[3] = threeBillion + 18;
        remoteGids[4] = threeBillion + 20;

        comm->barrier ();
        cerr << myRank << ": Calling getRemoteIndexList" << endl;
        comm->barrier ();
        map->getRemoteIndexList (remoteGids (), remotePids (), remoteLids ());

        TEST_COMPARE_ARRAYS( remotePids (), expectedRemotePids () );
        TEST_COMPARE_ARRAYS( remoteLids (), expectedRemoteLids () );
      }
      catch (std::exception& e) {
        locallyThrew = 1;
        std::ostringstream os;
        os << "Proc " << myRank << ": getRemoteIndexList locally threw "
          "exception: " << e.what () << endl;
        cerr << os.str ();
      }
    }
    else if (myRank == 1) {
      try {
        Array<int> expectedRemotePids (numRemoteGids);
        std::fill (expectedRemotePids.begin (), expectedRemotePids.end (), 0);
        Array<int> expectedRemoteLids (numRemoteGids);
        expectedRemoteLids[0] = 0;
        expectedRemoteLids[1] = 1;
        expectedRemoteLids[2] = 2;
        expectedRemoteLids[3] = 3;
        expectedRemoteLids[4] = 4;
        expectedRemoteLids[5] = 5;
        remoteGids[0] = globalFirstGid;
        remoteGids[1] = threeBillion;
        remoteGids[2] = threeBillion + 2;
        remoteGids[3] = threeBillion + 4;
        remoteGids[4] = threeBillion + 6;
        remoteGids[5] = threeBillion + 8;

        comm->barrier ();
        cerr << myRank << ": Calling getRemoteIndexList" << endl;
        comm->barrier ();
        map->getRemoteIndexList (remoteGids (), remotePids (), remoteLids ());

        TEST_COMPARE_ARRAYS( remotePids (), expectedRemotePids () );
        TEST_COMPARE_ARRAYS( remoteLids (), expectedRemoteLids () );
      }
      catch (std::exception& e) {
        locallyThrew = 1;
        std::ostringstream os;
        os << "Proc " << myRank << ": getRemoteIndexList locally threw "
          "exception: " << e.what () << endl;
        cerr << os.str ();
      }
    }

    reduceAll<int, int> (*comm, REDUCE_MAX, locallyThrew, outArg (globallyThrew));
    TEST_EQUALITY_CONST( globallyThrew, 0 );
  }

  cout << endl; // make TimeMonitor output neat on test line
  Teuchos::TimeMonitor::summarize();

  int globalSuccess_int = -1;
  Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
  TEST_EQUALITY_CONST( globalSuccess_int, 0 );
}
} // (anonymous)
