// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_TestingUtilities.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include "Tpetra_Map.hpp"
#include <algorithm> // std::min

namespace { // (anonymous)

using Teuchos::as;
using Teuchos::Comm;
using Teuchos::MpiComm;
using Teuchos::null;
using Teuchos::outArg;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::REDUCE_MIN;
using Teuchos::REDUCE_SUM;
using Teuchos::reduceAll;
using Teuchos::toString;
using std::endl;
typedef Tpetra::global_size_t GST;

// Tests that apply for any subset Map.  Must be called collectively
// over subsetComm (NOT the original Map's communicator!).
//
// NOTE: Don't call this with the unit test's original FancyOStream.
// That thing only prints to Process 0 in MPI_COMM_WORLD by default.
// Calls to this function may not necessarily include that process.
//
// NOTE: Don't call any functions (including Map methods on origMap)
// that may be collective over the original communicator!  This
// function is only meant to be called on processes that participate
// in the subset communicator (i.e., where subsetComm is nonnull).
// Calling collective functions on origMap or its communicator will
// likely result in a hang!  This includes Map methods like
// getRemoteIndexList and isOneToOne.
template<class LO, class GO, class NT>
void
testSubsetMapOverSubsetComm (int& gblSuccess, // output argument; 0 means false
                             std::ostream& out,
                             const Teuchos::RCP<const Tpetra::Map<LO, GO, NT> >& origMap,
                             const Teuchos::RCP<const Tpetra::Map<LO, GO, NT> >& subsetMap,
                             const Teuchos::RCP<const Teuchos::Comm<int> >& subsetComm)
{
  int lclSuccess = 1; // to be modified below

  // We need to establish Map validity collectively.  If the input Map
  // is not valid, we need to be able to bail out early, so that we
  // won't hang when making calls with collective semantics on the Map
  // itself.  Thus, we need the subset communicator to be nonnull, so
  // that we can call collectives on it.
  TEUCHOS_TEST_FOR_EXCEPTION
    (subsetComm.is_null (), std::logic_error, "Input subsetComm is null.  "
     "This should never happen.");

  // Rank of the calling process in the original communicator.  This
  // is always sensible, since the original Map's communicator
  // includes all processes in the subset Map.
  const int myRank = origMap->getComm ()->getRank ();

  if (subsetMap.is_null ()) {
    lclSuccess = 0;
    out << "Process " << myRank << ": Null Map" << endl;
  }
  else if (subsetMap->getComm ().is_null ()) {
    lclSuccess = 0;
    out << "Process " << myRank << ": Nonnull Map with null communicator"
        << endl;
  }
  if (subsetMap->getComm ().getRawPtr () != subsetComm.getRawPtr ()) {
    lclSuccess = 0;
    out << "Process " << myRank << ": Nonnull Map with nonnull communicator, "
      "but communicator is not the same (pointer) as the input subset "
      "communicator" << endl;
  }

  reduceAll<int, int> (*subsetComm, REDUCE_MIN,
                       lclSuccess, outArg (gblSuccess));
  if (gblSuccess != 1) {
    return; // subset Map is in no shape to continue
  }

  // At this point, the subset Map is sensible, so we can at least
  // test its features.  Let's start with the obviously local features
  // of the Map.  We'll all-reduce after checking those, just to make
  // sure that we get quick output in case something goes wrong
  // thereafter.  Then, we'll compare their local and global indices
  // on the calling process.  This may segfault or have a bounds
  // checking error if the global-to-local or local-to-global Maps
  // were not set up correctly, so it pays to do an all-reduce right
  // before, just so that we get _some_ output.  Finally, we'll do
  // another all-reduce, then check the features of the Map that may
  // require communication (and thus have collective semantics).

  if (subsetMap->getLocalNumElements () != origMap->getLocalNumElements ()) {
    lclSuccess = 0;
    out << "subsetMap->getLocalNumElements() = "
        << subsetMap->getLocalNumElements ()
        << " != origMap->getLocalNumElements() = "
        << origMap->getLocalNumElements ()
        << endl;
  }
  if (subsetMap->getMinLocalIndex () != origMap->getMinLocalIndex ()) {
    lclSuccess = 0;
    out << "subsetMap->getMinLocalIndex() = "
        << subsetMap->getMinLocalIndex ()
        << " != origMap->getMinLocalIndex() = "
        << origMap->getMinLocalIndex ()
        << endl;
  }
  if (subsetMap->getMaxLocalIndex () != origMap->getMaxLocalIndex ()) {
    lclSuccess = 0;
    out << "subsetMap->getMaxLocalIndex() = "
        << subsetMap->getMaxLocalIndex ()
        << " != origMap->getMaxLocalIndex() = "
        << origMap->getMaxLocalIndex ()
        << endl;
  }
  if (subsetMap->getMinGlobalIndex () != origMap->getMinGlobalIndex ()) {
    lclSuccess = 0;
    out << "subsetMap->getMinGlobalIndex() = "
        << subsetMap->getMinGlobalIndex ()
        << " != origMap->getMinGlobalIndex() = "
        << origMap->getMinGlobalIndex ()
        << endl;
  }
  if (subsetMap->getMaxGlobalIndex () != origMap->getMaxGlobalIndex ()) {
    lclSuccess = 0;
    out << "subsetMap->getMaxGlobalIndex() = "
        << subsetMap->getMaxGlobalIndex ()
        << " != origMap->getMaxGlobalIndex() = "
        << origMap->getMaxGlobalIndex ()
        << endl;
  }

  reduceAll<int, int> (*subsetComm, REDUCE_MIN,
                       lclSuccess, outArg (gblSuccess));
  if (gblSuccess != 1) {
    return; // we know the subset Map is messed up
  }

  try {
    const LO lclNumInds =
      std::min (static_cast<LO> (subsetMap->getLocalNumElements ()),
                static_cast<LO> (origMap->getLocalNumElements ()));
    Teuchos::Array<LO> badLclInds;
    Teuchos::Array<std::pair<GO, GO> > badGblInds;
    bool foundBadInd = false;
    for (LO lclInd = 0; lclInd < lclNumInds; ++lclInd) {
      const GO newGblInd = subsetMap->getGlobalElement (lclInd);
      const GO oldGblInd = origMap->getGlobalElement (lclInd);
      if (newGblInd != oldGblInd) {
        foundBadInd = true;
        badLclInds.push_back (lclInd);
        badGblInds.push_back (std::make_pair (oldGblInd, newGblInd));
      }
    }

    if (foundBadInd) {
      lclSuccess = 0;
      out << "The following local indices have global indices in the two Maps "
        "that do not match: " << Teuchos::toString (badLclInds)
          << endl << "Here are their corresponding global indices in "
        "(origMap, subsetMap): ";
      // out << Teuchos::toString (badGblInds) << endl; // doesn't compile, alas
      out << "[";
      const LO numBad = static_cast<LO> (badGblInds.size ());
      for (LO k = 0; k < numBad; ++k) {
        out << "(" << badGblInds[k].first << "," << badGblInds[k].second << ")";
        if (k + 1 != numBad) {
          out << ", ";
        }
      }
      out << "]" << endl;
    }
  }
  catch (std::exception& e) {
    lclSuccess = 0;
    out << "Process " << myRank << ": Testing local and global indices threw "
      "an exception: " << e.what () << endl;
  }
  catch (...) {
    lclSuccess = 0;
    out << "Process " << myRank << ": Testing local and global indices threw "
      "an exception, not a subclass of std::exception" << endl;
  }

  reduceAll<int, int> (*subsetComm, REDUCE_MIN,
                       lclSuccess, outArg (gblSuccess));
  if (gblSuccess != 1) {
    return; // we know the subset Map is messed up
  }

  // Use the original Map to compute how many entries the subset Map
  // _should_ have.  This is just the sum of all the local entry
  // counts in the original Map, over the subset comm.
  const GO lclNumEnt = static_cast<GO> (origMap->getLocalNumElements ());
  GO gblNumEnt = 0; // output argument
  reduceAll<int, GO> (*subsetComm, REDUCE_SUM, lclNumEnt, outArg (gblNumEnt));

  if (gblNumEnt != static_cast<GO> (subsetMap->getGlobalNumElements ())) {
    lclSuccess = 0;
    out << "subsetMap->getGlobalElements() = "
        << subsetMap->getGlobalNumElements ()
        << " != " << gblNumEnt << endl;
  }

  reduceAll<int, int> (*subsetComm, REDUCE_MIN,
                       lclSuccess, outArg (gblSuccess));
  // if (gblSuccess != 1) {
  //   return; // we know the subset Map is messed up
  // }
}

// Tests that apply for any subset Map.  Must be called collectively
// over the original Map's communicator.
//
// NOTE: Don't call this unless all the tests in
// testSubsetMapOverSubsetComm (see above) passed.
//
// NOTE: Don't call this with the unit test's original FancyOStream.
// That thing only prints to Process 0 in MPI_COMM_WORLD by default.
// We want to collect all the test's output.
template<class LO, class GO, class NT>
void
testSubsetMapOverOrigComm (int& gblSuccess, // output argument; 0 means false
                           std::ostream& out,
                           const Teuchos::RCP<const Tpetra::Map<LO, GO, NT> >& origMap,
                           const Teuchos::RCP<const Tpetra::Map<LO, GO, NT> >& subsetMap)
{
  int lclSuccess = 1; // to be modified below

  Teuchos::RCP<const Teuchos::Comm<int> > origComm = origMap->getComm ();
  Teuchos::RCP<const Teuchos::Comm<int> > subsetComm = subsetMap->getComm ();
  const int myRank = origComm->getRank ();

  // This is collective over the original comm, not the subset.
  const bool origMapIsOneToOne = origMap->isOneToOne ();
  bool subsetMapIsOneToOne = origMapIsOneToOne; // will revise below

  // isOneToOne is collective, so only call it on processes that
  // participate in the subset comm.
  if (! subsetComm.is_null ()) {
    // isOneToOne may do multiple rounds of communication, so even if
    // we catch exceptions, deadlock may still happen.  We just have
    // to do our best.
    try {
      subsetMapIsOneToOne = subsetMap->isOneToOne ();
    }
    catch (std::exception& e) {
      lclSuccess = 0;
      out << "Process " << myRank << ": subsetMap->isOneToOne() threw an "
        "exception: " << e.what () << endl;
    }
    catch (...) {
      lclSuccess = 0;
      out << "Process " << myRank << ": subsetMap->isOneToOne() threw an "
        "exception, not a subclass of std::exception" << endl;
    }
  }

  // If the original Map is one to one, then the subset Map must also
  // be one to one.  However, if the original Map is _not_ one to one,
  // the subset Map may or may not be, depending on what processes got
  // excluded.
  if (origMapIsOneToOne && ! subsetMapIsOneToOne) {
    lclSuccess = 0;
    out << "Original Map is one to one, but subset Map is NOT one to one"
        << endl;
  }

  reduceAll<int, int> (*origComm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  if (gblSuccess != 1) {
    return; // subset Map is in no shape to continue
  }
}

// Tests that apply for any subset Map.  Must be called collectively
// over the original Map's communicator.
//
// This function calls both testSubsetMapOverSubsetComm and
// testSubsetMapOverOrigComm (see above for both).
//
// NOTE: Don't call this with the unit test's original FancyOStream.
// That thing only prints to Process 0 in MPI_COMM_WORLD by default.
// We want to collect all the test's output.
template<class LO, class GO, class NT>
void
testSubsetMap (int& gblSuccess, // output argument; 0 means false
               std::ostream& out,
               const Teuchos::RCP<const Tpetra::Map<LO, GO, NT> >& origMap,
               const Teuchos::RCP<const Teuchos::Comm<int> >& origComm,
               const Teuchos::RCP<const Tpetra::Map<LO, GO, NT> >& subsetMap,
               const Teuchos::RCP<const Teuchos::Comm<int> >& subsetComm)
{
  int lclSuccess = 1; // to be modified below
  std::ostringstream errStrm; // for collecting test output

  if (! subsetComm.is_null ()) {
    testSubsetMapOverSubsetComm (gblSuccess, errStrm, origMap,
                                 subsetMap, subsetComm);
    lclSuccess = gblSuccess;
  }
  // The subset comm may not necessarily include Process 0 in the
  // original comm, so do (another) all-reduce over the original comm
  // in order to collect the first test's error state and output.
  reduceAll<int, int> (*origComm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  if (gblSuccess != 1) {
    Tpetra::Details::gathervPrint (out, errStrm.str (), *origComm);
    return; // we know that subset Map is messed up, so don't continue
  }

  testSubsetMapOverOrigComm (gblSuccess, errStrm, origMap, subsetMap);
  if (gblSuccess != 1) {
    Tpetra::Details::gathervPrint (out, errStrm.str (), *origComm);
    return;
  }
}

// This test is only meaningful in an MPI build.
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( Map, replaceCommWithSubset, LO, GO, NT )
{
  typedef Tpetra::Map<LO, GO, NT> map_type;
  typedef typename Teuchos::Array<GO>::size_type size_type;
  int lclSuccess = 1; // to be modified below
  int gblSuccess = 0; // output argument
  std::ostringstream errStrm; // for collecting output

  out << "Test Tpetra::Map::replaceCommWithSubset" << endl;
  Teuchos::OSTab tab1 (out);

  // We really want MPI_COMM_WORLD here, not just the test's default
  // communicator.  This will make sure that every process that _can_
  // participate, _does_ participate.
  out << "Create original Comm that wraps MPI_COMM_WORLD" << endl;
  RCP<const Comm<int> > origComm = rcp (new MpiComm<int> (MPI_COMM_WORLD));
  const int numProcs = origComm->getSize ();
  const int myRank = origComm->getRank ();

  out << "Create original Map, in which all processes have a nonzero number "
    "of entries" << endl;
  const size_type numGidsPerProc = 3;
  const size_type myNumGids = numGidsPerProc;
  Teuchos::Array<GO> myGids (myNumGids);
  for (size_type k = 0; k < myNumGids; ++k) {
    myGids[k] = as<GO> (myRank) *
      as<GO> (numGidsPerProc) +
      as<GO> (k);
  }
  const GST globalNumElts =
    static_cast<GST> (numGidsPerProc) * static_cast<GST> (numProcs);
  const GO indexBase = 0;
  RCP<const map_type> origMap (new map_type (globalNumElts, myGids (),
                                             indexBase, origComm));

  out << "Split the original communicator into {Process 0}, "
    "{all other processes}" << endl;
  // If color == 0 on a process, then that process will participate in
  // the subset communicator.  Otherwise, it won't.
  const int color = (myRank == 0) ? 0 : 1;
  const int key = 0;
  RCP<const Comm<int> > subsetComm;
  try {
    subsetComm = origComm->split (color, key);
  }
  catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Process " << myRank << ": origComm->split(...) threw an "
      "exception: " << e.what () << endl;
  }
  catch (...) {
    lclSuccess = 0;
    errStrm << "Process " << myRank << ": origComm->split(...) threw an "
      "exception not a subclass of std::exception" << endl;
  }

  // The new communicator should be nonnull on all processes.
  // Make sure of that before we continue.
  if (subsetComm.is_null ()) {
    lclSuccess = 0;
    errStrm << "Process " << myRank << ": result of splitting the original "
      "comm is null." << endl;
  }
  reduceAll<int, int> (*origComm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY( gblSuccess, 1 );
  if (gblSuccess != 1) {
    Tpetra::Details::gathervPrint (out, errStrm.str (), *origComm);
    return; // no sense in continuing
  }

  // replaceCommWithSubset must be called collectively on the
  // _original_ communicator.  It leaves the original Map and its
  // (original) communicator unchanged.
  out << "Call replaceCommWithSubset with the new communicator "
    "to create the new Map" << endl;
  RCP<const map_type> subsetMap;
  try {
    subsetMap = origMap->replaceCommWithSubset (subsetComm);
  }
  catch (std::exception& e) {
    lclSuccess = 0;
    errStrm << "Process " << myRank << ": "
      "origComm->replaceCommWithSubset(subsetComm) threw an exception: "
            << e.what () << endl;
  }
  catch (...) {
    lclSuccess = 0;
    errStrm << "Process " << myRank << ": "
      "origMap->replaceCommWithSubset(subsetComm) threw an exception "
      "not a subclass of std::exception" << endl;
  }
  reduceAll<int, int> (*origComm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY( gblSuccess, 1 );
  if (gblSuccess != 1) {
    Tpetra::Details::gathervPrint (out, errStrm.str (), *origComm);
    return; // no sense in continuing
  }

  out << "replaceCommWithSubset didn't throw.  Now test null/nonnull." << endl;

  if (subsetMap.is_null ()) {
    lclSuccess = 0;
    errStrm << "Process " << myRank << ": replaceCommWithSubset returned a "
      "null Map" << endl;
  }
  else if (subsetMap->getComm ().is_null ()) {
    lclSuccess = 0;
    errStrm << "Process " << myRank << ": replaceCommWithSubset returned a "
      "nonnull Map with a null communicator" << endl;
  }

  reduceAll<int, int> (*origComm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY( gblSuccess, 1 );
  if (gblSuccess != 1) {
    Tpetra::Details::gathervPrint (out, errStrm.str (), *origComm);
    return; // no sense in continuing
  }

  testSubsetMap (gblSuccess, out, origMap, origComm, subsetMap, subsetComm);

  // if (subsetMap->getMinAllGlobalIndex () != origMap->getMinGlobalIndex ()) {
  //   lclSuccess = 0;
  //   err << "subsetMap->getMinAllGlobalIndex() = "
  //       << subsetMap->getMinAllGlobalIndex ()
  //       << " != origMap->getMinGlobalIndex() = "
  //       << origMap->getMinGlobalIndex ()
  //       << endl;
  // }
  // if (subsetMap->getMaxAllGlobalIndex () != origMap->getMaxGlobalIndex ()) {
  //   lclSuccess = 0;
  //   err << "subsetMap->getMaxAllGlobalIndex() = "
  //       << subsetMap->getMaxAllGlobalIndex ()
  //       << " != origMap->getMaxGlobalIndex() = "
  //       << origMap->getMaxGlobalIndex ()
  //       << endl;
  // }

  TEST_EQUALITY_CONST( gblSuccess, 1 );
}

//
// INSTANTIATIONS OF TESTS
//

#define UNIT_TEST_GROUP( LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Map, replaceCommWithSubset, LO, GO, NODE )

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_LGN( UNIT_TEST_GROUP )

} // namespace (anonymous)



