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

#include "TpetraNew_Map.hpp"
#include "Tpetra_TestingUtilities.hpp"

// Exercise Map::getMyGlobalIndices().

namespace { // (anonymous)

using Tpetra::TestingUtilities::getDefaultComm;
using Teuchos::Comm;
using Teuchos::outArg;
using Teuchos::RCP;
using Teuchos::REDUCE_MIN;
using Teuchos::reduceAll;
using std::endl;

using map_type = TpetraNew::Map;
using LO = map_type::local_ordinal_type;
using GO = map_type::global_ordinal_type;

void
testGids (bool& success,
          Teuchos::FancyOStream& out,
          const map_type& gblMap)
{
  const GO gblInvalid = Tpetra::Details::OrdinalTraits<GO>::invalid ();

  try {
    auto gblInds = gblMap.getMyGlobalIndices ();
    auto gblInds_host = Kokkos::create_mirror_view (gblInds);
    Kokkos::deep_copy (gblInds_host, gblInds);

    TEST_EQUALITY( static_cast<size_t> (gblInds_host.size ()),
                   static_cast<size_t> (gblMap.getNodeNumElements ()) );
    if (success) {
      const LO numLclElts = static_cast<LO> (gblInds_host.size ());

      // Test the reported global indices.
      for (LO lid = 0; lid < numLclElts; ++lid) {
        const GO expectedGid = gblMap.getGlobalElement (lid);
        const GO reportedGid = gblInds_host(lid);

        // Make sure that the (global) Map behaves as expected.
        TEST_INEQUALITY( expectedGid, gblInvalid );
        // Make sure gblInds_host contains the right global index.
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
TEUCHOS_UNIT_TEST( getMyGlobalIndices, UniformContig )
{
  int lclSuccess = 1;
  int gblSuccess = 1;

  Teuchos::OSTab tab0 (out);
  out << "Map::getMyGlobalIndices: Uniform Contiguous Map" << endl;
  Teuchos::OSTab tab1 (out);

  RCP<const Comm<int> > comm = getDefaultComm ();
  const LO numLclElts = 10;
  const GO numGblElts = GO (comm->getSize ()) * GO (numLclElts);
  const GO indexBase = 0;

  map_type gblMap (numGblElts, indexBase, comm);

  // Sanity check on the global Map.
  TEST_EQUALITY( gblMap.getNodeNumElements (), LO (numLclElts) );
  TEST_EQUALITY( gblMap.getMinLocalIndex (), LO (0) );
  TEST_EQUALITY( gblMap.getMaxLocalIndex (), LO (numLclElts - 1) );

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
TEUCHOS_UNIT_TEST( getMyGlobalIndices, NonuniformContig )
{
  int lclSuccess = 1;
  int gblSuccess = 1;

  Teuchos::OSTab tab0 (out);
  out << "Map::getMyGlobalIndices: Nonuniform Contiguous Map" << endl;
  Teuchos::OSTab tab1 (out);

  RCP<const Comm<int> > comm = getDefaultComm ();
  const int myRank = comm->getRank ();

  // Process p gets 5+p indices.
  const LO numLclElts = static_cast<LO> (5 + myRank);
  const GO numGblElts = Tpetra::Details::OrdinalTraits<GO>::invalid ();
  const GO indexBase = 0;

  map_type gblMap (numGblElts, numLclElts, indexBase, comm);

  // Sanity check on the global Map.
  TEST_EQUALITY( gblMap.getNodeNumElements (), LO (numLclElts) );
  TEST_EQUALITY( gblMap.getMinLocalIndex (), LO (0) );
  TEST_EQUALITY( gblMap.getMaxLocalIndex (), LO (numLclElts - 1) );

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
TEUCHOS_UNIT_TEST( getMyGlobalIndices, Noncontig )
{
  int lclSuccess = 1;
  int gblSuccess = 1;

  Teuchos::OSTab tab0 (out);
  out << "Map::getMyGlobalIndices: Noncontiguous Map" << endl;
  Teuchos::OSTab tab1 (out);

  RCP<const Comm<int> > comm = getDefaultComm ();
  const int myRank = comm->getRank ();

  // Process p gets 5 indices.
  const LO numLclElts = 5;
  const GO numGblElts = Tpetra::Details::OrdinalTraits<GO>::invalid ();
  const GO indexBase = 0;

  Teuchos::Array<GO> eltList (numLclElts);
  const GO myGidStart = GO (GO (myRank + 1) * numLclElts - 1);
  for (LO lid = 0; lid < numLclElts; ++lid) {
    eltList[lid] = static_cast<GO> (myGidStart - lid);
  }

  map_type gblMap (numGblElts, eltList (), indexBase, comm);

  // Sanity check on the global Map.
  TEST_EQUALITY( gblMap.getNodeNumElements (), LO (numLclElts) );
  TEST_EQUALITY( gblMap.getMinLocalIndex (), LO (0) );
  TEST_EQUALITY( gblMap.getMaxLocalIndex (), LO (numLclElts - 1) );

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

} // namespace (anonymous)


