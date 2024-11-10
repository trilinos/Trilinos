// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Map.hpp"
#include "Tpetra_TestingUtilities.hpp"
#include <type_traits> // std::is_same

namespace {

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::reduceAll;
  using Teuchos::outArg;
  using Tpetra::global_size_t;
  using std::endl;

  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Map, getMinGlobalIndex_nonuniform, LO, GO )
  {
    out << "Test: Map, getMinGlobalIndex" << std::endl;
    Teuchos::OSTab tab0 (out);

    // create a comm
    auto comm = Tpetra::getDefaultComm();
    const int numRanks = comm->getSize();
    const int myRank = comm->getRank();
    // create a contiguous uniform distributed map with numLocal entries per process
    const size_t        numLocalRef = 5;
    const size_t        numLocal    = (myRank == 0) ? 0 : numLocalRef;
    const global_size_t numGlobal   = (numRanks - 1)*numLocalRef;
    const GO indexBase = 0;

    Tpetra::Map<LO,GO> map (numGlobal, numLocal, indexBase, comm);
    // create data to check validity of the map
    const GO myMinGID = (myRank == 0) ?
      std::numeric_limits<GO>::max () :
      static_cast<GO> ((myRank - 1)*numLocal);
    const GO myMaxGID = (myRank == 0) ?
      std::numeric_limits<GO>::lowest () :
      static_cast<GO> (myRank*numLocal - 1);
    GO minAllGIDs, maxAllGIDs;
    reduceAll<int, GO>(*comm, Teuchos::REDUCE_MIN, myMinGID, outArg(minAllGIDs));
    reduceAll<int, GO>(*comm, Teuchos::REDUCE_MAX, myMaxGID, outArg(maxAllGIDs));
    //
    TEST_EQUALITY(map.getGlobalNumElements(), numGlobal);
    TEST_EQUALITY(map.getLocalNumElements(), numLocal);
    TEST_EQUALITY(map.getIndexBase(), indexBase);
    TEST_EQUALITY(map.getMinGlobalIndex(), myMinGID);
    TEST_EQUALITY(map.getMaxGlobalIndex(), myMaxGID );
    TEST_EQUALITY(map.getMinAllGlobalIndex(), minAllGIDs);
    TEST_EQUALITY(map.getMaxAllGlobalIndex(), maxAllGIDs);
    // All procs fail if any proc fails
    int gblSuccess = -1;
    reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(gblSuccess) );
    TEST_EQUALITY_CONST( gblSuccess, 0 );
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Map, getMinGlobalIndex_noncontiguous, LO, GO )
  {
    out << "Test: Map, getMinGlobalIndex" << std::endl;
    Teuchos::OSTab tab0 (out);

    // create a comm
    auto comm = Tpetra::getDefaultComm();
    const int numRanks = comm->getSize();
    const int myRank = comm->getRank();
    // create a contiguous uniform distributed map with numLocal entries per node
    const size_t numLocalRef = 5;
    const size_t numLocal = (myRank == 0) ? static_cast<size_t> (0) : numLocalRef;
    const global_size_t numGlobal = static_cast<global_size_t>((numRanks - 1) * numLocalRef);
    const GO indexBase = 0;
    const GO actualBase = 1;
    //
    Teuchos::Array<GO> GIDs;
    if(myRank > 0) {
      GIDs.resize(numLocal);
      GIDs[0] = actualBase + myRank*numLocal;
      for (size_t i = 1; i < numLocal; ++i) {
        GIDs[i] = GIDs[i-1]+1;
      }
    }

    Tpetra::Map<LO,GO> map (numGlobal, GIDs (), indexBase, comm);
    // create data to check validity of the map
    const GO myMinGID = (myRank == 0) ?
      std::numeric_limits<GO>::max() :
      static_cast<GO>(actualBase + myRank*numLocal);
    const GO myMaxGID = (myRank == 0) ?
      std::numeric_limits<GO>::lowest() :
      static_cast<GO>(actualBase + (myRank + 1)*numLocal - 1);
    GO minAllGIDs, maxAllGIDs;
    reduceAll<int, GO>(*comm, Teuchos::REDUCE_MIN, myMinGID, outArg(minAllGIDs));
    reduceAll<int, GO>(*comm, Teuchos::REDUCE_MAX, myMaxGID, outArg(maxAllGIDs));
    //
    TEST_EQUALITY(map.getGlobalNumElements(), numGlobal);
    TEST_EQUALITY(map.getLocalNumElements(), numLocal);
    TEST_EQUALITY(map.getIndexBase(), indexBase);
    TEST_EQUALITY(map.getMinGlobalIndex(), myMinGID);
    TEST_EQUALITY(map.getMaxGlobalIndex(), myMaxGID );
    TEST_EQUALITY(map.getMinAllGlobalIndex(), minAllGIDs);
    TEST_EQUALITY(map.getMaxAllGlobalIndex(), maxAllGIDs);
    Teuchos::ArrayView<const GO> glist = map.getLocalElementList();
    TEST_COMPARE_ARRAYS( map.getLocalElementList(), GIDs);
    // All procs fail if any proc fails
    int gblSuccess = -1;
    reduceAll (*comm, Teuchos::REDUCE_SUM, success ? 0 : 1,
	       outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 0 );
  }

  //
  // INSTANTIATIONS
  //

#define UNIT_TEST_GROUP( LO, GO ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, getMinGlobalIndex_nonuniform, LO, GO ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Map, getMinGlobalIndex_noncontiguous, LO, GO )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_LG(UNIT_TEST_GROUP)

} // namespace (anonymous)
