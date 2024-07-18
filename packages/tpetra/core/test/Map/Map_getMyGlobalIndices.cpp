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

// Exercise Tpetra::Map::getMyGlobalIndices().

namespace { // (anonymous)

using Tpetra::TestingUtilities::getDefaultComm;
using Teuchos::Comm;
using Teuchos::outArg;
using Teuchos::RCP;
using Teuchos::REDUCE_MIN;
using Teuchos::reduceAll;
using std::endl;
typedef Tpetra::global_size_t GST;

template<class MapType>
void
testGids (bool& success,
          Teuchos::FancyOStream& out,
          const MapType& gblMap)
{
  typedef typename MapType::local_ordinal_type LO;
  typedef typename MapType::global_ordinal_type GO;
  const GO gblInvalid = Tpetra::Details::OrdinalTraits<GO>::invalid ();

  try {
    auto gblInds = gblMap.getMyGlobalIndices ();

    TEST_EQUALITY( static_cast<size_t> (gblInds.size ()),
                   static_cast<size_t> (gblMap.getLocalNumElements ()) );
    if (success) {
      const LO numLclElts = static_cast<LO> (gblInds.size ());

      // Test the reported global indices.
      for (LO lid = 0; lid < numLclElts; ++lid) {
        const GO expectedGid = gblMap.getGlobalElement (lid);
        const GO reportedGid = gblInds(lid);

        // Make sure that the (global) Map behaves as expected.
        TEST_INEQUALITY( expectedGid, gblInvalid );
        // Make sure gblInds contains the right global index.
        TEST_EQUALITY( expectedGid, reportedGid );
      }
    }
  }
  catch (...) {
    success = false;
  }
  // No need for reduction here; this function's caller does that.
}

//
// Test with a uniform contiguous Map.
//
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( getMyGlobalIndices, UniformContig, LO, GO, NT )
{
  typedef Tpetra::Map<LO, GO, NT> map_type;
  int lclSuccess = 1;
  int gblSuccess = 1;

  Teuchos::OSTab tab0 (out);
  out << "Map::getMyGlobalIndices: Uniform Contiguous Map" << endl;
  Teuchos::OSTab tab1 (out);

  RCP<const Comm<int> > comm = getDefaultComm ();
  const LO numLclElts = 10;

  const GST numGblElts = static_cast<GST> (comm->getSize () * numLclElts);
  const GO indexBase = 0;

  map_type gblMap (numGblElts, indexBase, comm);

  // Sanity check on the global Map.
  TEST_EQUALITY( gblMap.getLocalNumElements (), static_cast<size_t> (numLclElts) );
  TEST_EQUALITY( gblMap.getMinLocalIndex (), static_cast<LO> (0) );
  TEST_EQUALITY( gblMap.getMaxLocalIndex (), static_cast<LO> (numLclElts - 1) );

  testGids (success, out, gblMap);

  // Make sure that all processes succeeded.
  lclSuccess = success ? 1 : 0;
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY_CONST( gblSuccess, 1 );

  if (gblSuccess != 1) {
    out << "The test failed on at least one process." << endl;
    // Make sure that the test fails, even if Process 0 was fine.
    success = false;
  }
}

//
// Test with a NONuniform contiguous Map.
//
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( getMyGlobalIndices, NonuniformContig, LO, GO, NT )
{
  typedef Tpetra::Map<LO, GO, NT> map_type;
  int lclSuccess = 1;
  int gblSuccess = 1;

  Teuchos::OSTab tab0 (out);
  out << "Map::getMyGlobalIndices: Nonuniform Contiguous Map" << endl;
  Teuchos::OSTab tab1 (out);

  RCP<const Comm<int> > comm = getDefaultComm ();
  const int myRank = comm->getRank ();

  // Process p gets 5+p indices.
  const LO numLclElts = static_cast<LO> (5 + myRank);
  const GST numGblElts = Tpetra::Details::OrdinalTraits<GST>::invalid ();
  const GO indexBase = 0;

  map_type gblMap (numGblElts, static_cast<size_t> (numLclElts), indexBase, comm);

  // Sanity check on the global Map.
  TEST_EQUALITY( gblMap.getLocalNumElements (), static_cast<size_t> (numLclElts) );
  TEST_EQUALITY( gblMap.getMinLocalIndex (), static_cast<LO> (0) );
  TEST_EQUALITY( gblMap.getMaxLocalIndex (), static_cast<LO> (numLclElts - 1) );

  testGids (success, out, gblMap);

  // Make sure that all processes succeeded.
  lclSuccess = success ? 1 : 0;
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY_CONST( gblSuccess, 1 );

  if (gblSuccess != 1) {
    out << "The test failed on at least one process." << endl;
    // Make sure that the test fails, even if Process 0 was fine.
    success = false;
    return;
  }
}

//
// Test with a NONcontiguous Map.
//
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( getMyGlobalIndices, Noncontig, LO, GO, NT )
{
  typedef Tpetra::Map<LO, GO, NT> map_type;
  int lclSuccess = 1;
  int gblSuccess = 1;

  Teuchos::OSTab tab0 (out);
  out << "Map::getMyGlobalIndices: Noncontiguous Map" << endl;
  Teuchos::OSTab tab1 (out);

  RCP<const Comm<int> > comm = getDefaultComm ();
  const int myRank = comm->getRank ();

  // Process p gets 5 indices.
  const LO numLclElts = static_cast<LO> (5);
  const GST numGblElts = Tpetra::Details::OrdinalTraits<GST>::invalid ();
  const GO indexBase = 0;

  Teuchos::Array<GO> eltList (numLclElts);
  const GO myGidStart = static_cast<GO> ((myRank + 1) * numLclElts  - 1);
  for (LO lid = 0; lid < numLclElts; ++lid) {
    eltList[lid] = static_cast<GO> (myGidStart - lid);
  }

  map_type gblMap (numGblElts, eltList (), indexBase, comm);

  // Sanity check on the global Map.
  TEST_EQUALITY( gblMap.getLocalNumElements (), static_cast<size_t> (numLclElts) );
  TEST_EQUALITY( gblMap.getMinLocalIndex (), static_cast<LO> (0) );
  TEST_EQUALITY( gblMap.getMaxLocalIndex (), static_cast<LO> (numLclElts - 1) );

  testGids (success, out, gblMap);

  // Make sure that all processes succeeded.
  lclSuccess = success ? 1 : 0;
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY_CONST( gblSuccess, 1 );

  if (gblSuccess != 1) {
    out << "The test failed on at least one process." << endl;
    // Make sure that the test fails, even if Process 0 was fine.
    success = false;
  }
}

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( getMyGlobalIndices, UniformContig, LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( getMyGlobalIndices, NonuniformContig, LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( getMyGlobalIndices, Noncontig, LO, GO, NT )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_LGN(UNIT_TEST_GROUP)

} // namespace (anonymous)


