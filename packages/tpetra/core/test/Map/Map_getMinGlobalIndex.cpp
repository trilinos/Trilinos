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
#include <Tpetra_Map.hpp>
#include <Tpetra_TestingUtilities.hpp>
#include <type_traits> // std::is_same

// FINISH: add testing of operator==, operator!=, operator=, copy construct
// put these into test_same_as and test_is_compatible

namespace {

  using Teuchos::null;

  using Tpetra::TestingUtilities::getNode;
  using Tpetra::TestingUtilities::getDefaultComm;

  using Tpetra::createUniformContigMapWithNode;

  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::arcp;
  using Teuchos::rcp;
  using Teuchos::outArg;
  using Teuchos::Tuple;
  using Teuchos::tuple;
  using Tpetra::Map;
  using Tpetra::global_size_t;
  using Tpetra::DefaultPlatform;
  using std::sort;
  using std::find;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::Comm;
  using std::endl;

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Map, getMinGlobalIndex_nonuniform, LO, GO )
  {
    out << "Test: Map, getMinGlobalIndex" << std::endl;
    Teuchos::OSTab tab0 (out);

    // create a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numRanks = comm->getSize();
    const int myRank = comm->getRank();
    // create a contiguous uniform distributed map with numLocal entries per node
    const size_t        numLocalRef = 5;
    const size_t        numLocal    = (myRank == 0) ? 0 : numLocalRef;
    const global_size_t numGlobal   = (numRanks - 1)*numLocalRef;
    const GO indexBase = 0;
    //

    Map<LO,GO> map (numGlobal, numLocal, indexBase, comm);
    // create data to check validity of the map
    const GO myMinGID = (myRank == 0) ?
      std::numeric_limits<GO>::max () :
      static_cast<GO> ((myRank - 1)*numLocal);
    const GO myMaxGID = (myRank == 0) ?
      std::numeric_limits<GO>::lowest () :
      static_cast<GO> (myRank*numLocal - 1);
    GO minAllGIDs, maxAllGIDs;
    Teuchos::reduceAll<int, GO>(*comm, Teuchos::REDUCE_MIN, myMinGID, outArg(minAllGIDs));
    Teuchos::reduceAll<int, GO>(*comm, Teuchos::REDUCE_MAX, myMaxGID, outArg(maxAllGIDs));
    //
    TEST_EQUALITY(map.getGlobalNumElements(), numGlobal);
    TEST_EQUALITY(map.getNodeNumElements(), numLocal);
    TEST_EQUALITY(map.getIndexBase(), indexBase);
    TEST_EQUALITY(map.getMinGlobalIndex(), myMinGID);
    TEST_EQUALITY(map.getMaxGlobalIndex(), myMaxGID );
    TEST_EQUALITY(map.getMinAllGlobalIndex(), minAllGIDs);
    TEST_EQUALITY(map.getMaxAllGlobalIndex(), maxAllGIDs);
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Map, getMinGlobalIndex_noncontiguous, LO, GO )
  {
    out << "Test: Map, getMinGlobalIndex" << std::endl;
    Teuchos::OSTab tab0 (out);

    // create a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    const int numRanks = comm->getSize();
    const int myRank = comm->getRank();
    // create a contiguous uniform distributed map with numLocal entries per node
    const size_t numLocalRef = 5;
    const size_t numLocal = (myRank == 0) ? static_cast<size_t> (0) : numLocalRef;
    const global_size_t numGlobal = static_cast<global_size_t>((numRanks - 1) * numLocalRef);
    const GO indexBase = 0;
    const GO actualBase = 1;
    //
    Array<GO> GIDs;
    if(myRank > 0) {
      GIDs.resize(numLocal);
      GIDs[0] = actualBase + myRank*numLocal;
      for (size_t i = 1; i < numLocal; ++i) {
        GIDs[i] = GIDs[i-1]+1;
      }
    }

    Map<LO,GO> map (numGlobal, GIDs (), indexBase, comm);
    // create data to check validity of the map
    const GO myMinGID = (myRank == 0) ?
      std::numeric_limits<GO>::max() :
      static_cast<GO>(actualBase + myRank*numLocal);
    const GO myMaxGID = (myRank == 0) ?
      std::numeric_limits<GO>::lowest() :
      static_cast<GO>(actualBase + (myRank + 1)*numLocal - 1);
    GO minAllGIDs, maxAllGIDs;
    Teuchos::reduceAll<int, GO>(*comm, Teuchos::REDUCE_MIN, myMinGID, outArg(minAllGIDs));
    Teuchos::reduceAll<int, GO>(*comm, Teuchos::REDUCE_MAX, myMaxGID, outArg(maxAllGIDs));
    //
    TEST_EQUALITY(map.getGlobalNumElements(), numGlobal);
    TEST_EQUALITY(map.getNodeNumElements(), numLocal);
    TEST_EQUALITY(map.getIndexBase(), indexBase);
    TEST_EQUALITY(map.getMinGlobalIndex(), myMinGID);
    TEST_EQUALITY(map.getMaxGlobalIndex(), myMaxGID );
    TEST_EQUALITY(map.getMinAllGlobalIndex(), minAllGIDs);
    TEST_EQUALITY(map.getMaxAllGlobalIndex(), maxAllGIDs);
    ArrayView<const GO> glist = map.getNodeElementList();
    TEST_COMPARE_ARRAYS( map.getNodeElementList(), GIDs);
    // All procs fail if any proc fails
    int globalSuccess_int = -1;
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, success ? 0 : 1, outArg(globalSuccess_int) );
    TEST_EQUALITY_CONST( globalSuccess_int, 0 );
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
