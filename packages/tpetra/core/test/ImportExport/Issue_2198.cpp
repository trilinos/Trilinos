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

//
// Test for Github Issue #2198 (feature request for new Import constructor).
//

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include "Teuchos_DefaultComm.hpp"

using Teuchos::outArg;
using Teuchos::rcp;
using Teuchos::REDUCE_MIN;
using Teuchos::reduceAll;
using std::endl;
typedef Tpetra::global_size_t GST;

template<class LO, class GO, class NT>
bool
importsLocallySame (Teuchos::FancyOStream& out,
                    const Tpetra::Import<LO, GO, NT>& X,
                    const Tpetra::Import<LO, GO, NT>& Y)
{
  Teuchos::OSTab tab0 (out);

  bool someSourceMapIsNull = false;
  bool sourceMapsSame = true;
  if (X.getSourceMap ().is_null ()) {
    someSourceMapIsNull = true;
    if (! Y.getSourceMap ().is_null ()) {
      sourceMapsSame = false;
    }
  }
  else if (Y.getSourceMap ().is_null ()) {
    someSourceMapIsNull = true;
    if (! X.getSourceMap ().is_null ()) {
      sourceMapsSame = false;
    }
  }
  else if (! X.getSourceMap ()->isSameAs (* (Y.getSourceMap ()))) {
    sourceMapsSame = false;
  }
  else if (! Y.getSourceMap ()->isSameAs (* (X.getSourceMap ()))) {
    sourceMapsSame = false;
  }

  if (someSourceMapIsNull) {
    out << "At least one source Map is null!" << endl;
  }
  if (! sourceMapsSame) {
    out << "Source Maps not same!" << endl;
    return false;
  }

  bool someTargetMapIsNull = false;
  bool targetMapsSame = true;
  if (X.getTargetMap ().is_null ()) {
    someTargetMapIsNull = true;
    if (! Y.getTargetMap ().is_null ()) {
      targetMapsSame = false;
    }
  }
  else if (Y.getTargetMap ().is_null ()) {
    someTargetMapIsNull = true;
    if (! X.getTargetMap ().is_null ()) {
      targetMapsSame = false;
    }
  }
  else if (! X.getTargetMap ()->isSameAs (* (Y.getTargetMap ()))) {
    targetMapsSame = false;
  }
  else if (! Y.getTargetMap ()->isSameAs (* (X.getTargetMap ()))) {
    targetMapsSame = false;
  }

  if (someTargetMapIsNull) {
    out << "At least one target Map is null!" << endl;
  }
  if (! targetMapsSame) {
    out << "Target Maps not same!" << endl;
    return false;
  }

  if (X.getNumSameIDs () != Y.getNumSameIDs ()) {
    out << "X.getNumSameIDs() = " << X.getNumSameIDs ()
        << " != Y.getNumSameIDs() = " << Y.getNumSameIDs () << endl;
    return false;
  }
  if (X.getNumPermuteIDs () != Y.getNumPermuteIDs ()) {
    out << "X.getNumPermuteIDs() = " << X.getNumPermuteIDs ()
        << " != Y.getNumPermuteIDs() = " << Y.getNumPermuteIDs () << endl;
    return false;
  }
  if (X.getNumRemoteIDs () != Y.getNumRemoteIDs ()) {
    out << "X.getNumRemoteIDs() = " << X.getNumRemoteIDs ()
        << " != Y.getNumRemoteIDs() = " << Y.getNumRemoteIDs () << endl;
    return false;
  }
  if (X.getNumExportIDs () != Y.getNumExportIDs ()) {
    out << "X.getNumExportIDs() = " << X.getNumExportIDs ()
        << " != Y.getNumExportIDs() = " << Y.getNumExportIDs () << endl;
    return false;
  }

  if (! std::equal (X.getPermuteFromLIDs (), Y.getPermuteFromLIDs ())) {
    out << "X.getPermuteFromLIDs() = " << Teuchos::toString (X.getPermuteFromLIDs ())
        << " != Y.getPermuteFromLIDs() = " << Teuchos::toString (Y.getPermuteFromLIDs ())
        << endl;
    return false;
  }
  if (! std::equal (X.getPermuteToLIDs (), Y.getPermuteToLIDs ())) {
    out << "X.getPermuteToLIDs() = " << Teuchos::toString (X.getPermuteToLIDs ())
        << " != Y.getPermuteToLIDs() = " << Teuchos::toString (Y.getPermuteToLIDs ())
        << endl;
    return false;
  }

  if (! std::equal (X.getExportLIDs (), Y.getExportLIDs ())) {
    out << "X.getExportLIDs() = " << Teuchos::toString (X.getExportLIDs ())
        << " != Y.getExportLIDs() = " << Teuchos::toString (Y.getExportLIDs ())
        << endl;
    return false;
  }
  if (! std::equal (X.getRemoteLIDs (), Y.getRemoteLIDs ())) {
    out << "X.getRemoteLIDs() = " << Teuchos::toString (X.getRemoteLIDs ())
        << " != Y.getRemoteLIDs() = " << Teuchos::toString (Y.getRemoteLIDs ())
        << endl;
    return false;
  }
  if (! std::equal (X.getRemotePIDs (), Y.getRemotePIDs ())) {
    out << "X.getRemotePIDs() = " << Teuchos::toString (X.getRemotePIDs ())
        << " != Y.getRemotePIDs() = " << Teuchos::toString (Y.getRemotePIDs ())
        << endl;
    return false;
  }

  return true;
}

template<class LO, class GO, class NT>
struct Issue2198TestInput {
  Teuchos::RCP<const Teuchos::Comm<int> > comm;
  std::vector<GO> sourceMapGlobalIndices;
  std::vector<GO> targetMapGlobalIndices;
  std::vector<int> targetMapProcessRanks;
  const GST globalNumSourceMapGlobalIndices;
  const GO indexBase;
};

template<class LO, class GO, class NT>
Teuchos::RCP<const Tpetra::Map<LO, GO, NT> >
makeSourceMapFromTestInput (const Issue2198TestInput<LO, GO, NT>& test)
{
  typedef Tpetra::Map<LO, GO, NT> map_type;
  return rcp (new map_type (test.globalNumSourceMapGlobalIndices,
                            test.sourceMapGlobalIndices.data (),
                            test.sourceMapGlobalIndices.size (),
                            test.indexBase,
                            test.comm));
}

template<class LO, class GO, class NT>
Teuchos::RCP<const Tpetra::Map<LO, GO, NT> >
makeUnoptimizedTargetMapFromTestInput (const Issue2198TestInput<LO, GO, NT>& test)
{
  typedef Tpetra::Map<LO, GO, NT> map_type;
  return rcp (new map_type (Teuchos::OrdinalTraits<GST>::invalid (),
                            test.targetMapGlobalIndices.data (),
                            test.targetMapGlobalIndices.size (),
                            test.indexBase,
                            test.comm));
}

template<class LO, class GO, class NT>
Issue2198TestInput<LO, GO, NT>
makeTest_A (const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
{
  const int myRank = comm->getRank ();

  std::vector<GO> sourceMapGlobalIndices;
  if (myRank == 0) {
    sourceMapGlobalIndices = {0, 4, 6};
  }
  else if (myRank == 1) {
    sourceMapGlobalIndices = {1, 2, 5};
  }
  else if (myRank == 2) {
    sourceMapGlobalIndices = {11, 14, 15};
  }
  else if (myRank == 3) {
    sourceMapGlobalIndices = {10, 12, 16};
  }
  else if (myRank == 4) {
    sourceMapGlobalIndices = {3, 7, 8, 9, 13};
  }

  std::vector<GO> targetMapGlobalIndices;
  if (myRank == 0) {
    targetMapGlobalIndices = {1, 3, 7};
  }
  else if (myRank == 1) {
    targetMapGlobalIndices = {3, 9, 10};
  }
  else if (myRank == 2) {
    targetMapGlobalIndices = {6, 7, 13};
  }
  else if (myRank == 3) {
    targetMapGlobalIndices = {9, 13, 15};
  }
  else if (myRank == 4) {
    targetMapGlobalIndices = {4, 5, 11, 12};
  }

  std::vector<int> targetMapProcessRanks;
  if (myRank == 0) {
    targetMapProcessRanks = {1, 4, 4};
  }
  else if (myRank == 1) {
    targetMapProcessRanks = {4, 4, 3};
  }
  else if (myRank == 2) {
    targetMapProcessRanks = {0, 4, 4};
  }
  else if (myRank == 3) {
    targetMapProcessRanks = {4, 4, 2};
  }
  else if (myRank == 4) {
    targetMapProcessRanks = {0, 1, 2, 3};
  }

  const GST globalNumSourceMapGlobalIndices = 17;
  const GO indexBase = 0;
  return {comm, sourceMapGlobalIndices,
      targetMapGlobalIndices,
      targetMapProcessRanks,
      globalNumSourceMapGlobalIndices,
      indexBase};
}


TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( ImportExport, Issue2198_A, GO, NT )
{
  typedef Tpetra::Map<>::local_ordinal_type LO;
  //typedef Tpetra::Map<LO, GO, NT> map_type;
  typedef Tpetra::Import<LO, GO, NT> import_type;
  int lclSuccess = 1;
  int gblSuccess = 1;

  out << "Test for Issue 2198 (new Import constructor)" << endl;
  Teuchos::OSTab tab1 (out);

  Teuchos::RCP<const Teuchos::Comm<int> > origComm = Teuchos::DefaultComm<int>::getComm ();
  const int origNumProcs = origComm->getSize ();
  constexpr int minNumProcs = 5;
  constexpr int neededNumProcs = 5;
  if (origNumProcs < minNumProcs) {
    out << "Test FAILED; must be run on at least " << minNumProcs << " MPI "
      "processes, but the test's input communicator has size " << origNumProcs
      << endl;
    success = false;
    return;
  }
  const int origMyRank = origComm->getRank ();
  const int color = (origMyRank < neededNumProcs) ? 0 : 1;
  auto comm = origComm->split (color, origMyRank);
  if (color == 1) {
    return;
  }
  //const int myRank = comm->getRank ();
  // Only processes in comm -- that is, in the first 0
  // .. neededNumProcs-1 processes of origComm -- participate in what
  // follows.

  Issue2198TestInput<LO, GO, NT> test_A = makeTest_A<LO, GO, NT> (comm);
  auto sourceMap = makeSourceMapFromTestInput (test_A);
  auto expectedUnoptimizedTargetMap = makeUnoptimizedTargetMapFromTestInput (test_A);

  {
    constexpr bool optimized = false;
    import_type expectedUnoptimizedImport (sourceMap, expectedUnoptimizedTargetMap);
    import_type actualImport (sourceMap,
                              test_A.targetMapGlobalIndices.data (),
                              test_A.targetMapProcessRank.data (),
                              test_A.targetMapGlobalIndices.size (),
                              optimized);
    const bool impsLocallySame =
      importsLocallySame (out, expectedUnoptimizedImport, actualImport);
    TEST_ASSERT( impsLocallySame );
  }

  lclSuccess = success ? 1 : 0;
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY_CONST( gblSuccess, 1 );
}

