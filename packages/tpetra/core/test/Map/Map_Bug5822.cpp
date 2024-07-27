// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_Map.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_UnitTestHarness.hpp>

namespace { // (anonymous)
using Tpetra::global_size_t;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::as;
using Teuchos::Comm;
using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::tuple;
using std::endl;

// Test requires exactly 2 MPI processes
TEUCHOS_UNIT_TEST( Map, Bug5822_StartWith3Billion )
{
  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm ();
  const int numProcs = comm->getSize ();
  const int myRank = comm->getRank ();
  TEUCHOS_TEST_FOR_EXCEPTION(numProcs != 2, std::logic_error,
    "This test only makes sense to run with 2 MPI processes.");

#ifdef HAVE_TPETRA_INT_LONG_LONG
  typedef long long GO;
  if (sizeof (long long) <= 4) {
    out << "sizeof (long long) = " << sizeof (long long) << " <= 4.  "
      "This test only makes sense if sizeof (long long) >= 8, "
      "since the test is supposed to exercise GIDs > 2 billion.";
    return;
  }
#else // NOT HAVE_TPETRA_INT_LONG_LONG
  typedef long GO;
  if (sizeof (long) <= 4) {
    out << "sizeof (long) = " << sizeof (long) << " <= 4.  "
      "This test only makes sense if sizeof (long) >= 8, "
      "since the test is supposed to exercise GIDs > 2 billion.";
    return;
  }
#endif // HAVE_TPETRA_INT_LONG_LONG
  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Details::DefaultTypes::node_type NT;
  typedef Tpetra::Map<LO, GO, NT> map_type;

  const size_t localNumElts = 5;
  const global_size_t globalNumElts = comm->getSize () * localNumElts;
  const GO globalFirstGid = 3000000000L;
  const GO localFirstGid = globalFirstGid + 2 * localNumElts * comm->getRank ();

  // Make sure that the GIDs are not contiguous, in case Map tries to
  // detect contiguous GIDs for us.  The point is to exercise a
  // noncontiguous Map with GIDs all greater than 2 billion (so that
  // they would overflow a 32-bit integer).
  Array<GO> myGids (localNumElts);
  Array<GO> myGidsExpected (localNumElts);
  for (size_t k = 0; k < localNumElts; ++k) {
    myGids[k] = localFirstGid + as<GO> (2*k);
    // Make a copy, just to make sure that Map's constructor didn't
    // sneakily change the input ArrayView.
    myGidsExpected[k] = myGids[k];
  }
  // Proc 0: myGids = [3B, 3B+2, 3B+4, 3B+6, 3B+8].
  // Proc 1: myGids = [3B+10, 3B+12, 3B+14, 3B+16, 3B+18].

  // Tpetra::Map requires that the index base and the first GID be equal.
  const GO indexBase = globalFirstGid;
  RCP<const map_type> map (new map_type (globalNumElts, myGids (), indexBase, comm));

  ArrayView<const GO> myGidsFound = map->getLocalElementList ();
  TEST_COMPARE_ARRAYS( myGidsExpected (), myGidsFound () );

  Array<GO> remoteGids (5);
  Array<int> remotePids (5, -1);
  Array<LO> remoteLids (5, Teuchos::OrdinalTraits<LO>::invalid ());
  if (myRank == 0) {
    Array<int> expectedRemotePids (5);
    std::fill (expectedRemotePids.begin (), expectedRemotePids.end (), 1);
    Array<int> expectedRemoteLids (5);
    expectedRemoteLids[0] = 0;
    expectedRemoteLids[1] = 1;
    expectedRemoteLids[2] = 2;
    expectedRemoteLids[3] = 3;
    expectedRemoteLids[4] = 4;
    remoteGids[0] = globalFirstGid + 10;
    remoteGids[1] = globalFirstGid + 12;
    remoteGids[2] = globalFirstGid + 14;
    remoteGids[3] = globalFirstGid + 16;
    remoteGids[4] = globalFirstGid + 18;
    map->getRemoteIndexList (remoteGids (), remotePids (), remoteLids ());

    TEST_COMPARE_ARRAYS( remotePids (), expectedRemotePids () );
    TEST_COMPARE_ARRAYS( remoteLids (), expectedRemoteLids () );
  } else if (myRank == 1) {
    Array<int> expectedRemotePids (5);
    std::fill (expectedRemotePids.begin (), expectedRemotePids.end (), 0);
    Array<int> expectedRemoteLids (5);
    expectedRemoteLids[0] = 0;
    expectedRemoteLids[1] = 1;
    expectedRemoteLids[2] = 2;
    expectedRemoteLids[3] = 3;
    expectedRemoteLids[4] = 4;
    remoteGids[0] = globalFirstGid;
    remoteGids[1] = globalFirstGid + 2;
    remoteGids[2] = globalFirstGid + 4;
    remoteGids[3] = globalFirstGid + 6;
    remoteGids[4] = globalFirstGid + 8;
    map->getRemoteIndexList (remoteGids (), remotePids (), remoteLids ());

    TEST_COMPARE_ARRAYS( remotePids (), expectedRemotePids () );
    TEST_COMPARE_ARRAYS( remoteLids (), expectedRemoteLids () );
  }
}
} // (anonymous)
