// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Tpetra_TestingUtilities.hpp>
#include <Tpetra_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_DefaultSerialComm.hpp>
#include <Tpetra_Map.hpp>

#ifdef HAVE_TPETRA_EXPLICIT_INSTANTIATION
#include <Tpetra_Map_def.hpp>
#include <Tpetra_Directory_def.hpp>
#endif // HAVE_TPETRA_EXPLICIT_INSTANTIATION

namespace { // (anonymous)

using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::as;
using Teuchos::Comm;
using Teuchos::MpiComm;
using Teuchos::outArg;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::REDUCE_MIN;
using Teuchos::reduceAll;
using Teuchos::SerialComm;
using Teuchos::toString;
using Teuchos::tuple;
using std::cerr;
using std::endl;

typedef Tpetra::global_size_t GST;

// This test is only meaningful in an MPI build.
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Map, removeEmptyProcesses_MpiComm_noncontigMap, LO, GO, NT)
{
  typedef Tpetra::Map<LO, GO, NT> map_type;
  typedef typename Array<GO>::size_type size_type;

  RCP<const Comm<int> > origComm = rcp(new MpiComm<int>(MPI_COMM_WORLD));
  const int numProcs = origComm->getSize();
  const int myRank = origComm->getRank();

  // All processes but the last one(myRank == numProcs-1) contain a
  // nonzero number of elements.
  const size_type numGidsPerProc = 3; // not counting the last process
  const size_type myNumGids = (myRank == numProcs - 1) ? 0 : numGidsPerProc;
  Array<GO> myGids(myNumGids);
  if (myNumGids > 0) {
    for (size_type k = 0; k < myNumGids; ++k) {
      myGids[k] = as<GO>(myRank) * as<GO>(numGidsPerProc) + as<GO>(k);
    }
  }

  const GST globalNumElts = as<GST>(numGidsPerProc) * as<GST>(numProcs - 1);
  const GO indexBase = 0;
  RCP<const map_type> origMap(new map_type(globalNumElts, myGids(),
                                           indexBase, origComm));
  RCP<const map_type> newMap = origMap->removeEmptyProcesses();

  // Test collectively for success, so the test doesn't hang on failure.
  int localSuccess = 1;
  std::ostringstream err;
  if (myRank == numProcs - 1) {
    if (! newMap.is_null()) {
      localSuccess = 0;
      err << "removeEmptyProcesses() should have returned null, but did not."
          << endl;
    }
  } else {
    if (newMap.is_null()) {
      localSuccess = 0;
      err << "removeEmptyProcesses() should not have returned null, but did."
          << endl;
    } else {
      RCP<const Comm<int> > newComm = newMap->getComm();
      if (newComm->getSize() != numProcs - 1) {
        localSuccess = 0;
        err << "New communicator should have " << (numProcs - 1)
            << " processes, but has " << newComm->getSize()
            << " processes instead." << endl;
      }

      if (newMap->getGlobalNumElements() != origMap->getGlobalNumElements()) {
        localSuccess = 0;
        err << "New Map has " << newMap->getGlobalNumElements() << " global "
            << "elements, which differs from the original number "
            << origMap->getGlobalNumElements() << "." << endl;
      }

      if (newMap->getLocalNumElements() != origMap->getLocalNumElements()) {
        localSuccess = 0;
        err << "New Map has " << newMap->getLocalNumElements() << " local "
          "elements, which differs from the original number "
          << origMap->getLocalNumElements() << "." << endl;
      }

      if (newMap->getIndexBase() != origMap->getIndexBase()) {
        localSuccess = 0;
        err << "New Map has index base " << newMap->getIndexBase()
            << ", which differs from the original Map's index base "
          << origMap->getIndexBase() << "." << endl;
      }

      if (newMap->getMinLocalIndex() != origMap->getMinLocalIndex()) {
        localSuccess = 0;
        err << "New Map has min local index " << newMap->getMinLocalIndex()
            << ", which differs from the original Map's min local index "
            << origMap->getMinLocalIndex() << "." << endl;
      }

      if (newMap->getMaxLocalIndex() != origMap->getMaxLocalIndex()) {
        localSuccess = 0;
        err << "New Map has max local index " << newMap->getMaxLocalIndex()
            << ", which differs from the original Map's max local index "
            << origMap->getMaxLocalIndex() << "." << endl;
      }

      if (newMap->getMinGlobalIndex() != origMap->getMinGlobalIndex()) {
        localSuccess = 0;
        err << "New Map has local min global index "
            << newMap->getMinGlobalIndex()
            << ", which differs from the original Map's local min global index "
            << origMap->getMinGlobalIndex() << "." << endl;
      }

      if (newMap->getMaxGlobalIndex() != origMap->getMaxGlobalIndex()) {
        localSuccess = 0;
        err << "New Map has local max global index "
            << newMap->getMaxGlobalIndex()
            << ", which differs from the original Map's local max global index "
            << origMap->getMaxGlobalIndex() << "." << endl;
      }

      if (newMap->getMinAllGlobalIndex() != origMap->getMinAllGlobalIndex()) {
        localSuccess = 0;
        err << "New Map has global min global index "
            << newMap->getMinAllGlobalIndex() << ", which differs from the "
            << "original Map's global min global index "
            << origMap->getMinAllGlobalIndex() << "." << endl;
      }

      if (newMap->getMaxAllGlobalIndex() != origMap->getMaxAllGlobalIndex()) {
        localSuccess = 0;
        err << "New Map has global max global index "
            << newMap->getMaxAllGlobalIndex() << ", which differs from the "
            << "original Map's global max global index "
            << origMap->getMaxAllGlobalIndex() << "." << endl;
      }

      ArrayView<const GO> myNewGids = newMap->getLocalElementList();
      if (myNewGids.size() != myGids.size() ||
          ! std::equal(myNewGids.begin(), myNewGids.end(), myGids.begin())) {
        localSuccess = 0;
        err << "New Map has local GID list " << toString(myNewGids) << ", but "
            << "should have local GID list " << toString(myGids()) << "."
            << endl;
      }
    }
  }

  int globalSuccess = 0;
  reduceAll(*origComm, REDUCE_MIN, localSuccess, outArg(globalSuccess));
  if (globalSuccess == 0) {
    if (myRank == 0) {
      cerr << "TEST FAILED" << endl
           << "Error messages from each process:" << endl << endl;
    }
    for (int p = 0; p < numProcs; ++p) {
      if (myRank == p) {
        cerr << "Process " << myRank << ": " << err.str() << endl;
      }
      origComm->barrier(); // Give time for output to finish.
      origComm->barrier();
      origComm->barrier();
    }
  }
  TEST_EQUALITY(globalSuccess, 1);

  if (globalSuccess != 1){
    success = false;
    return;
  }
}

// This test is only meaningful in an MPI build.
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Map, removeEmptyProcesses_MpiComm_contigMap, LO, GO, NT)
{
  typedef Tpetra::Map<LO, GO, NT> map_type;

  RCP<const Comm<int> > origComm = rcp(new MpiComm<int>(MPI_COMM_WORLD));
  const int numProcs = origComm->getSize();
  const int myRank = origComm->getRank();

  // Let Map compute the global number of elements.
  const GST globalNumElts = Teuchos::OrdinalTraits<GST>::invalid();
  // All processes but the last one (myRank == numProcs-1) contain a
  // nonzero number of elements.
  const size_t numEltsPerProc = 3; // on all but the last process
  const size_t myNumElts = (myRank == numProcs - 1) ? 0 : numEltsPerProc;
  const GO indexBase = 0;
  RCP<const map_type> origMap(new map_type(globalNumElts, myNumElts,
                                             indexBase, origComm));
  RCP<const map_type> newMap = origMap->removeEmptyProcesses();

  // Test collectively for success, so the test doesn't hang on failure.
  int localSuccess = 1;
  std::ostringstream err;
  if (myRank == numProcs - 1) {
    if (! newMap.is_null()) {
      localSuccess = 0;
      err << "removeEmptyProcesses() should have returned null, but did not."
          << endl;
    }
  } else {
    if (newMap.is_null()) {
      localSuccess = 0;
      err << "removeEmptyProcesses() should not have returned null, but did."
          << endl;
    } else {
      RCP<const Comm<int> > newComm = newMap->getComm();
      if (newComm->getSize() != numProcs - 1) {
        localSuccess = 0;
        err << "New communicator should have " << (numProcs - 1)
            << " processes, but has " << newComm->getSize()
            << " processes instead." << endl;
      }

      if (newMap->getGlobalNumElements() != origMap->getGlobalNumElements()) {
        localSuccess = 0;
        err << "New Map has " << newMap->getGlobalNumElements() << " global "
            << "elements, which differs from the original number "
            << origMap->getGlobalNumElements() << "." << endl;
      }

      if (newMap->getLocalNumElements() != origMap->getLocalNumElements()) {
        localSuccess = 0;
        err << "New Map has " << newMap->getLocalNumElements() << " local "
          "elements, which differs from the original number "
          << origMap->getLocalNumElements() << "." << endl;
      }

      if (newMap->getIndexBase() != origMap->getIndexBase()) {
        localSuccess = 0;
        err << "New Map has index base " << newMap->getIndexBase()
            << ", which differs from the original Map's index base "
          << origMap->getIndexBase() << "." << endl;
      }

      if (newMap->getMinLocalIndex() != origMap->getMinLocalIndex()) {
        localSuccess = 0;
        err << "New Map has min local index " << newMap->getMinLocalIndex()
            << ", which differs from the original Map's min local index "
            << origMap->getMinLocalIndex() << "." << endl;
      }

      if (newMap->getMaxLocalIndex() != origMap->getMaxLocalIndex()) {
        localSuccess = 0;
        err << "New Map has max local index " << newMap->getMaxLocalIndex()
            << ", which differs from the original Map's max local index "
            << origMap->getMaxLocalIndex() << "." << endl;
      }

      if (newMap->getMinGlobalIndex() != origMap->getMinGlobalIndex()) {
        localSuccess = 0;
        err << "New Map has local min global index "
            << newMap->getMinGlobalIndex()
            << ", which differs from the original Map's local min global index "
            << origMap->getMinGlobalIndex() << "." << endl;
      }

      if (newMap->getMaxGlobalIndex() != origMap->getMaxGlobalIndex()) {
        localSuccess = 0;
        err << "New Map has local max global index "
            << newMap->getMaxGlobalIndex()
            << ", which differs from the original Map's local max global index "
            << origMap->getMaxGlobalIndex() << "." << endl;
      }

      if(newMap->getMinAllGlobalIndex() != origMap->getMinAllGlobalIndex()) {
        localSuccess = 0;
        err << "New Map has global min global index "
            << newMap->getMinAllGlobalIndex() << ", which differs from the "
            << "original Map's global min global index "
            << origMap->getMinAllGlobalIndex() << "." << endl;
      }

      if (newMap->getMaxAllGlobalIndex() != origMap->getMaxAllGlobalIndex()) {
        localSuccess = 0;
        err << "New Map has global max global index "
            << newMap->getMaxAllGlobalIndex() << ", which differs from the "
            << "original Map's global max global index "
            << origMap->getMaxAllGlobalIndex() << "." << endl;
      }

      ArrayView<const GO> myNewGids =
        newMap->getLocalElementList();
      ArrayView<const GO> myGids =
        origMap->getLocalElementList();
      if (myNewGids.size() != myGids.size() ||
          ! std::equal(myNewGids.begin(), myNewGids.end(), myGids.begin())) {
        localSuccess = 0;
        err << "New Map has local GID list " << toString(myNewGids) << ", but "
            << "should have local GID list " << toString(myGids()) << "."
            << endl;
      }
    }
  }

  int globalSuccess = 0;
  reduceAll(*origComm, REDUCE_MIN, localSuccess, outArg(globalSuccess));
  if (globalSuccess == 0) {
    if (myRank == 0) {
      cerr << "TEST FAILED" << endl
           << "Error messages from each process:" << endl << endl;
    }
    for (int p = 0; p < numProcs; ++p) {
      if (myRank == p) {
        cerr << "Process " << myRank << ": " << err.str() << endl;
      }
      origComm->barrier(); // Give time for output to finish.
      origComm->barrier();
      origComm->barrier();
    }
  }
  TEST_EQUALITY(globalSuccess, 1);

  if (globalSuccess != 1){
    success = false;
    return;
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( Map, removeEmptyProcesses_SerialComm1, LO, GO, NT )
{
  typedef Tpetra::Map<LO, GO, NT> map_type;
  typedef typename Array<GO>::size_type size_type;

  RCP<const Comm<int> > origComm = rcp(new SerialComm<int>);
  const int numProcs = origComm->getSize();
  const int myRank = origComm->getRank();

  // The calling process (there is only one in the communicator) contains a nonzero number of elements.
  const size_type numGidsPerProc = 3;
  const size_type myNumGids = numGidsPerProc;
  Array<GO> myGids(myNumGids);
  for (size_type k = 0; k < myNumGids; ++k) {
    myGids[k] = as<GO>(myRank) * as<GO>(numGidsPerProc) + as<GO>(k);
  }

  const GST globalNumElts = as<GST>(numGidsPerProc) * as<GST>(numProcs);
  const GO indexBase = 0;
  RCP<const map_type> origMap(new map_type(globalNumElts, myGids(),
                                             indexBase, origComm));
  RCP<const map_type> newMap = origMap->removeEmptyProcesses();

  // Test collectively for success, so the test doesn't hang on failure.
  int localSuccess = 1;
  std::ostringstream err;

  if (newMap.is_null()) {
    localSuccess = 0;
    err << "removeEmptyProcesses() should not have returned null, but did."
        << endl;
  } else {
    RCP<const Comm<int> > newComm = newMap->getComm();
    // No processes should be excluded in this case.
    if (newComm->getSize() != numProcs) {
      localSuccess = 0;
      err << "New communicator should have " << numProcs
          << " processes, but has " << newComm->getSize()
          << " processes instead." << endl;
    }

    if (newMap->getGlobalNumElements() != origMap->getGlobalNumElements()) {
      localSuccess = 0;
      err << "New Map has " << newMap->getGlobalNumElements() << " global "
          << "elements, which differs from the original number "
          << origMap->getGlobalNumElements() << "." << endl;
    }

    if (newMap->getLocalNumElements() != origMap->getLocalNumElements()) {
      localSuccess = 0;
      err << "New Map has " << newMap->getLocalNumElements() << " local "
        "elements, which differs from the original number "
          << origMap->getLocalNumElements() << "." << endl;
    }

    if (newMap->getIndexBase() != origMap->getIndexBase()) {
      localSuccess = 0;
      err << "New Map has index base " << newMap->getIndexBase()
          << ", which differs from the original Map's index base "
          << origMap->getIndexBase() << "." << endl;
    }

    if (newMap->getMinLocalIndex() != origMap->getMinLocalIndex()) {
      localSuccess = 0;
      err << "New Map has min local index " << newMap->getMinLocalIndex()
          << ", which differs from the original Map's min local index "
          << origMap->getMinLocalIndex() << "." << endl;
    }

    if (newMap->getMaxLocalIndex() != origMap->getMaxLocalIndex()) {
      localSuccess = 0;
      err << "New Map has max local index " << newMap->getMaxLocalIndex()
          << ", which differs from the original Map's max local index "
          << origMap->getMaxLocalIndex() << "." << endl;
    }

    if (newMap->getMinGlobalIndex() != origMap->getMinGlobalIndex()) {
      localSuccess = 0;
      err << "New Map has local min global index "
          << newMap->getMinGlobalIndex()
          << ", which differs from the original Map's local min global index "
          << origMap->getMinGlobalIndex() << "." << endl;
    }

    if (newMap->getMaxGlobalIndex() != origMap->getMaxGlobalIndex()) {
      localSuccess = 0;
      err << "New Map has local max global index "
          << newMap->getMaxGlobalIndex()
          << ", which differs from the original Map's local max global index "
          << origMap->getMaxGlobalIndex() << "." << endl;
    }

    if (newMap->getMinAllGlobalIndex() != origMap->getMinAllGlobalIndex()) {
      localSuccess = 0;
      err << "New Map has global min global index "
          << newMap->getMinAllGlobalIndex() << ", which differs from the "
          << "original Map's global min global index "
          << origMap->getMinAllGlobalIndex() << "." << endl;
    }

    if (newMap->getMaxAllGlobalIndex() != origMap->getMaxAllGlobalIndex()) {
      localSuccess = 0;
      err << "New Map has global max global index "
          << newMap->getMaxAllGlobalIndex() << ", which differs from the "
          << "original Map's global max global index "
          << origMap->getMaxAllGlobalIndex() << "." << endl;
    }

    ArrayView<const GO> myNewGids =
      newMap->getLocalElementList();
    if (myNewGids.size() != myGids.size() ||
        ! std::equal(myNewGids.begin(), myNewGids.end(), myGids.begin())) {
      localSuccess = 0;
      err << "New Map has local GID list " << toString(myNewGids) << ", but "
          << "should have local GID list " << toString(myGids()) << "."
          << endl;
    }
  }

  int globalSuccess = 0;
  reduceAll(*origComm, REDUCE_MIN, localSuccess, outArg(globalSuccess));
  if (globalSuccess == 0) {
    if (myRank == 0) {
      cerr << "TEST FAILED" << endl
           << "Error messages from each process:" << endl << endl;
    }
    for (int p = 0; p < numProcs; ++p) {
      if (myRank == p) {
        cerr << "Process " << myRank << ": " << err.str() << endl;
      }
      origComm->barrier(); // Give time for output to finish.
      origComm->barrier();
      origComm->barrier();
    }
  }
  TEST_EQUALITY(globalSuccess, 1);

  if (globalSuccess != 1){
    success = false;
    return;
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( Map, removeEmptyProcesses_SerialComm2, LO, GO, NT )
{
  typedef Tpetra::Map<LO, GO, NT> map_type;
  typedef typename Array<GO>::size_type size_type;

  RCP<const Comm<int> > origComm = rcp(new SerialComm<int>);
  const int numProcs = origComm->getSize();
  const int myRank = origComm->getRank();

  // The calling process (there is only one in the communicator) contains zero elements.
  const size_type numGidsPerProc = 0;
  const size_type myNumGids = numGidsPerProc;
  Array<GO> myGids(myNumGids);

  const GST globalNumElts = as<GST>(numGidsPerProc) *
    as<GST>(numProcs);
  const GO indexBase = 0;
  RCP<const map_type> origMap(new map_type(globalNumElts, myGids(),
                                             indexBase, origComm));
  RCP<const map_type> newMap = origMap->removeEmptyProcesses();

  // Test collectively for success, so the test doesn't hang on failure.
  int localSuccess = 1;
  std::ostringstream err;

  if (! newMap.is_null()) {
    localSuccess = 0;
    err << "removeEmptyProcesses() should have returned null, but did not."
        << endl;
  }

  int globalSuccess = 0;
  reduceAll(*origComm, REDUCE_MIN, localSuccess, outArg(globalSuccess));
  if (globalSuccess == 0) {
    if (myRank == 0) {
      cerr << "TEST FAILED" << endl
           << "Error messages from each process:" << endl << endl;
    }
    for (int p = 0; p < numProcs; ++p) {
      if (myRank == p) {
        cerr << "Process " << myRank << ": " << err.str() << endl;
      }
      origComm->barrier(); // Give time for output to finish.
      origComm->barrier();
      origComm->barrier();
    }
  }
  TEST_EQUALITY(globalSuccess, 1);

  if (globalSuccess != 1){
    success = false;
    return;
  }
}

//
// Instantiations of tests
//
#define UNIT_TEST_GROUP( LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Map, removeEmptyProcesses_MpiComm_noncontigMap, LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Map, removeEmptyProcesses_MpiComm_contigMap, LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Map, removeEmptyProcesses_SerialComm1, LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Map, removeEmptyProcesses_SerialComm2, LO, GO, NT )

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_LGN( UNIT_TEST_GROUP )

} // namespace (anonymous)
