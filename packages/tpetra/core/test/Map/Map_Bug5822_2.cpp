/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
// @HEADER
*/

// Some Macro Magic to ensure that if CUDA and KokkosCompat is enabled
// only the .cu version of this file is actually compiled
#include <Tpetra_ConfigDefs.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Tpetra_Map.hpp>
#include <Teuchos_TimeMonitor.hpp>

#ifdef HAVE_TPETRA_EXPLICIT_INSTANTIATION
#  include "Tpetra_Map_def.hpp"
#  include "Tpetra_Directory_def.hpp"
#  ifdef HAVE_TPETRA_FIXED_HASH_TABLE
#    include "Tpetra_Details_FixedHashTable_def.hpp"
#  else
#    include "Tpetra_HashTable_def.hpp"
#  endif // HAVE_TPETRA_FIXED_HASH_TABLE
#endif

using Tpetra::global_size_t;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::as;
using Teuchos::Comm;
using Teuchos::outArg;
using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::REDUCE_MAX;
using Teuchos::reduceAll;
using Teuchos::tuple;
using std::cerr;
using std::endl;
using std::cout;
using std::cin;

// This test works (and exercises the interesting case) in serial mode
// or for 1 MPI process, but it was originally written for 2 MPI
// processes.
TEUCHOS_UNIT_TEST( Map, Bug5822_StartWithZeroThenSkipTo3Billion )
{
  int locallyThrew = 0;
  int globallyThrew = 0;

  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm ();
  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();
  TEUCHOS_TEST_FOR_EXCEPTION(numProcs != 2, std::logic_error,
    "This test only makes sense to run with 2 MPI processes.");

#ifdef HAVE_TEUCHOS_LONG_LONG_INT
  typedef long long GO;
  if (sizeof (long long) <= 4) {
    out << "sizeof (long long) = " << sizeof (long long) << " <= 4.  "
      "This test only makes sense if sizeof (long long) >= 8, "
      "since the test is supposed to exercise GIDs > 2 billion.";
    return;
  }
#else // NOT HAVE_TEUCHOS_LONG_LONG_INT
  typedef long GO;
  if (sizeof (long) <= 4) {
    out << "sizeof (long) = " << sizeof (long) << " <= 4.  "
      "This test only makes sense if sizeof (long) >= 8, "
      "since the test is supposed to exercise GIDs > 2 billion.";
    return;
  }
#endif // HAVE_TEUCHOS_LONG_LONG_INT
  typedef int LO;
  typedef KokkosClassic::DefaultNode::DefaultNodeType NT;
  typedef Tpetra::Map<LO, GO, NT> map_type;

  // Proc 0 gets [0, 3B, 3B+2, 3B+4, 3B+8, 3B+10] (6 GIDs).
  // Proc 1 gets [3B+12, 3B+14, 3B+16, 3B+18, 3B+20] (5 GIDs).
  //
  // The GIDs are not contiguous in order to prevent Map from
  // detecting contiguity and bypassing the noncontiguous Map
  // case.
  const size_t localNumElts = (myRank == 0) ? 6 : 5;
  const global_size_t globalNumElts = 1 + 5*comm->getSize ();
  const GO globalFirstGid = 1L;
  const GO threeBillion = 3000000000L;

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
    myGids[0] = threeBillion + as<GO> ((localNumElts+1) * 2);
    myGidsExpected[0] = myGids[0];
    for (size_t k = 1; k < localNumElts; ++k) {
      myGids[k] = myGids[k-1] + 2;
      myGidsExpected[k] = myGids[k];
    }
  }

  RCP<NT> node;
  {
    Teuchos::ParameterList defaultParams;
    node = Teuchos::rcp (new NT (defaultParams));
  }
  // if (myRank == 0) {
  //   std::cout << "type '0' and hit enter" << std::endl;
  //   int zero;
  //   cin >> zero;
  // }
  // comm->barrier();
  // Tpetra::Map requires that the index base is less than the minimum GID.
  const GO indexBase = 0L;
  out << "Constructing the Map" << endl;
  RCP<const map_type> map;
  try {
    TEUCHOS_FUNC_TIME_MONITOR("Construct Map");
    map = Teuchos::rcp (new map_type (globalNumElts, myGids (), indexBase, comm, node));
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
    ArrayView<const GO> myGidsFound = map->getNodeElementList ();
    TEST_COMPARE_ARRAYS( myGidsExpected (), myGidsFound () );
  }

  cerr << myRank << ": Querying the Map for remote elements" << endl;
  // Proc 0 gets [0, 3B, 3B+2, 3B+4, 3B+8, 3B+10] (6 GIDs).
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
        remoteGids[0] = 0;
        remoteGids[1] = threeBillion;
        remoteGids[2] = threeBillion + 2;
        remoteGids[3] = threeBillion + 4;
        remoteGids[4] = threeBillion + 8;
        remoteGids[5] = threeBillion + 10;

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
}


