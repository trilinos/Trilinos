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
#include "Tpetra_Details_makeOptimizedColMap.hpp"
#include "Teuchos_DefaultComm.hpp"
#include <memory> // std::unique_ptr

using Teuchos::outArg;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::REDUCE_MIN;
using Teuchos::reduceAll;
using std::endl;
typedef Tpetra::global_size_t GST;

bool
trueEverywhere (const bool localTruthValue,
                const Teuchos::Comm<int>& comm)
{
  const int lcl = localTruthValue ? 1 : 0;
  int gbl = 0; // output argument
  Teuchos::reduceAll<int, int> (comm, Teuchos::REDUCE_MIN, lcl,
                                Teuchos::outArg (gbl));
  return gbl == 1;
}

// bool
// falseAnywhere (const bool localTruthValue,
//                const Teuchos::Comm<int>& comm)
// {
//   return ! trueEverywhere (localTruthValue, comm);
// }

template<class LO, class GO, class NT>
void
printMapCompactly (Teuchos::FancyOStream& out,
                   const Tpetra::Map<LO, GO, NT>& map)
{
  const Teuchos::Comm<int>& comm = * (map.getComm ());
  const int myRank = comm.getRank ();

  std::ostringstream lclOut;
  lclOut << "Proc " << myRank << ": [";
  if (map.getNodeNumElements () != 0) {
    for (LO lid = map.getMinLocalIndex (); lid <= map.getMaxLocalIndex (); ++lid) {
      const GO gid = map.getGlobalElement (lid);
      lclOut << gid;
      if (lid < map.getMaxLocalIndex ()) {
        lclOut << ", ";
      }
    }
  }
  lclOut << "]" << std::endl;
  Tpetra::Details::gathervPrint (out, lclOut.str (), comm);
}

template<class LO, class GO, class NT>
bool
importsLocallySame (Teuchos::FancyOStream& out,
                    const Tpetra::Import<LO, GO, NT>& X,
                    const std::string& X_name,
                    const Tpetra::Import<LO, GO, NT>& Y,
                    const std::string& Y_name,
                    const int myRank)
{
  Teuchos::OSTab tab0 (out);
  const std::string prefix = [&] () {
    std::ostringstream os;
    os << "Proc " << myRank << ": ";
    return os.str ();
  } ();

  bool same = true;
  if (X.getNumSameIDs () != Y.getNumSameIDs ()) {
    out << prefix << X_name << ".getNumSameIDs() = " << X.getNumSameIDs ()
        << " != " << Y_name << ".getNumSameIDs() = " << Y.getNumSameIDs () << endl;
    same = false;
  }
  // else {
  //   out << prefix << "Yay, numSames are the same!" << endl;
  // }

  if (X.getNumPermuteIDs () != Y.getNumPermuteIDs ()) {
    out << prefix << X_name << ".getNumPermuteIDs() = " << X.getNumPermuteIDs ()
        << " != " << Y_name << ".getNumPermuteIDs() = " << Y.getNumPermuteIDs () << endl;
    same = false;
  }
  // else {
  //   out << prefix << "Yay, numPermutes are the same!" << endl;
  // }

  if (X.getNumRemoteIDs () != Y.getNumRemoteIDs ()) {
    out << prefix << X_name << ".getNumRemoteIDs() = " << X.getNumRemoteIDs ()
        << " != " << Y_name << ".getNumRemoteIDs() = " << Y.getNumRemoteIDs () << endl;
    same = false;
  }
  // else {
  //   out << prefix << "Yay, numRemotes are the same!" << endl;
  // }

  if (X.getNumExportIDs () != Y.getNumExportIDs ()) {
    out << prefix << X_name << ".getNumExportIDs() = " << X.getNumExportIDs ()
        << " != " << Y_name << ".getNumExportIDs() = " << Y.getNumExportIDs () << endl;
    same = false;
  }
  // else {
  //   out << prefix << "Yay, numExports are the same!" << endl;
  // }

  if (! same) {
    return false; // can't do std::equal on sequences of different lengths
  }

  if (! std::equal (X.getPermuteFromLIDs ().begin (),
                    X.getPermuteFromLIDs ().end (),
                    Y.getPermuteFromLIDs ().begin ())) {
    out << prefix << X_name << ".getPermuteFromLIDs() = " << Teuchos::toString (X.getPermuteFromLIDs ())
        << " != " << Y_name << ".getPermuteFromLIDs() = " << Teuchos::toString (Y.getPermuteFromLIDs ())
        << endl;
    same = false;
  }
  // else {
  //   out << prefix << "Yay, permuteFroms are the same!" << endl;
  // }

  if (! std::equal (X.getPermuteToLIDs ().begin (),
                    X.getPermuteToLIDs ().end (),
                    Y.getPermuteToLIDs ().begin ())) {
    out << prefix << X_name << ".getPermuteToLIDs() = " << Teuchos::toString (X.getPermuteToLIDs ())
        << " != " << Y_name << ".getPermuteToLIDs() = " << Teuchos::toString (Y.getPermuteToLIDs ())
        << endl;
    same = false;
  }
  // else {
  //   out << prefix << "Yay, permuteTos are the same!" << endl;
  // }

  if (! std::equal (X.getExportLIDs ().begin (),
                    X.getExportLIDs ().end (),
                    Y.getExportLIDs ().begin ())) {
    out << prefix << X_name << ".getExportLIDs() = " << Teuchos::toString (X.getExportLIDs ())
        << " != " << Y_name << ".getExportLIDs() = " << Teuchos::toString (Y.getExportLIDs ())
        << endl;
    same = false;
  }
  // else {
  //   out << prefix << "Yay, exportLIDs are the same!" << endl;
  // }

  if (! std::equal (X.getRemoteLIDs ().begin (),
                    X.getRemoteLIDs ().end (),
                    Y.getRemoteLIDs ().begin ())) {
    out << prefix << X_name << ".getRemoteLIDs() = " << Teuchos::toString (X.getRemoteLIDs ())
        << " != " << Y_name << ".getRemoteLIDs() = " << Teuchos::toString (Y.getRemoteLIDs ())
        << endl;
    same = false;
  }
  // else {
  //   out << prefix << "Yay, remoteLIDs are the same!" << endl;
  // }

  return same;
}

template<class LO, class GO, class NT>
bool
importsGloballySame (Teuchos::FancyOStream& out,
                     const Tpetra::Import<LO, GO, NT>& X,
                     const std::string& X_name,
                     const Tpetra::Import<LO, GO, NT>& Y,
                     const std::string& Y_name)
{
  std::ostringstream lclOut;
  auto lclOutFOS = Teuchos::getFancyOStream (Teuchos::rcpFromRef (lclOut));
  const Teuchos::Comm<int>& comm = * (X.getSourceMap ()->getComm ());
  const int myRank = comm.getRank ();

  bool sourceMapsSame = true;
  if (! X.getSourceMap ()->isSameAs (* (Y.getSourceMap ()))) {
    sourceMapsSame = false;
  }
  else if (! Y.getSourceMap ()->isSameAs (* (X.getSourceMap ()))) {
    sourceMapsSame = false;
  }

  if (! sourceMapsSame) {
    out << "Source Maps not same!" << endl;
    return false;
  }
  else {
    out << "Yay, source Maps are the same!" << endl;
  }

  bool targetMapsSame = true;
  if (! X.getTargetMap ()->isSameAs (* (Y.getTargetMap ()))) {
    targetMapsSame = false;
  }
  else if (! Y.getTargetMap ()->isSameAs (* (X.getTargetMap ()))) {
    targetMapsSame = false;
  }

  if (! targetMapsSame) {
    out << "Target Maps not same!" << endl;
    return false;
  }
  else {
    out << "Yay, target Maps are the same!" << endl;
  }

  const bool lclImportsSame = importsLocallySame (*lclOutFOS, X, X_name, Y, Y_name, myRank);
  const bool gblImportsSame = trueEverywhere (lclImportsSame, comm);
  if (! gblImportsSame) {
    Tpetra::Details::gathervPrint (out, lclOut.str (), comm);
  }
  else {
    out << "Imports are the same on all processes!" << endl;
  }
  return gblImportsSame;
}

template<class LO, class GO, class NT>
struct Issue2198TestInput {
  Teuchos::RCP<const Teuchos::Comm<int> > comm;
  bool contiguous;
  LO localNumSourceMapGlobalIndices; // only valid if contiguous is true
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
  if (test.contiguous) {
    return rcp (new map_type (test.globalNumSourceMapGlobalIndices,
                              test.localNumSourceMapGlobalIndices,
                              test.indexBase,
                              test.comm));
  }
  else {
    return rcp (new map_type (test.globalNumSourceMapGlobalIndices,
                              test.sourceMapGlobalIndices.data (),
                              test.sourceMapGlobalIndices.size (),
                              test.indexBase,
                              test.comm));
  }
}

template<class LO, class GO, class NT>
Teuchos::RCP<const Tpetra::Map<LO, GO, NT> >
makeUnoptimizedTargetMapFromTestInput (const Issue2198TestInput<LO, GO, NT>& testInput,
                                       const ::Tpetra::Map<LO, GO, NT>& sourceMap)
{
  typedef Tpetra::Map<LO, GO, NT> map_type;

  const LO numLclSrcGids = sourceMap.getNodeNumElements ();
  const LO numInputGids = static_cast<LO> (testInput.targetMapGlobalIndices.size ());
  const LO numLclTgtGids = numLclSrcGids + numInputGids;
  std::vector<GO> tgtGids (numLclTgtGids);

  if (sourceMap.isContiguous ()) {
    GO curTgtGid = sourceMap.getMinGlobalIndex ();
    for (LO k = 0; k < numLclSrcGids; ++k, ++curTgtGid) {
      tgtGids[k] = curTgtGid;
    }
  }
  else {
    auto srcGids = sourceMap.getNodeElementList ();
    for (LO k = 0; k < numLclSrcGids; ++k) {
      tgtGids[k] = srcGids[k];
    }
  }
  std::copy (testInput.targetMapGlobalIndices.begin (),
             testInput.targetMapGlobalIndices.end (),
             tgtGids.begin () + numLclSrcGids);
  return rcp (new map_type (Teuchos::OrdinalTraits<GST>::invalid (),
                            tgtGids.data (),
                            tgtGids.size (),
                            testInput.indexBase,
                            testInput.comm));
}

template<class LO, class GO, class NT>
Issue2198TestInput<LO, GO, NT>
makeTest_A (const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
{
  const int myRank = comm->getRank ();
  constexpr bool contiguous = false;
  LO localNumSourceMapGlobalIndices = 0; // not used in this case
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
  return {comm, contiguous, localNumSourceMapGlobalIndices,
      sourceMapGlobalIndices,
      targetMapGlobalIndices,
      targetMapProcessRanks,
      globalNumSourceMapGlobalIndices,
      indexBase};
}

template<class LO, class GO, class NT>
Issue2198TestInput<LO, GO, NT>
makeTest_B (const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
{
  // Source Map: [0 1 2] [3 4 5 6] [7 8 9 10 11]
  //
  // Input for target Map: [5 11 3 8] [10 7] []
  //
  // Desired target Map, no sorting: [0 1 2, 5 11 3 8] [3 4 5 6, 10 7] [7 8 9 10 11]

  const int myRank = comm->getRank ();
  constexpr bool contiguous = true;
  LO localNumSourceMapGlobalIndices = 0;
  std::vector<GO> sourceMapGlobalIndices; // not used in this case
  if (myRank == 0) {
    localNumSourceMapGlobalIndices = 3;
  }
  else if (myRank == 1) {
    localNumSourceMapGlobalIndices = 4;
  }
  else if (myRank == 2) {
    localNumSourceMapGlobalIndices = 5;
  }

  std::vector<GO> targetMapGlobalIndices;
  if (myRank == 0) {
    targetMapGlobalIndices = {5, 11, 3, 8};
  }
  else if (myRank == 1) {
    targetMapGlobalIndices = {10, 7};
  }
  // else if (myRank == 2) {
  // // do nothing; empty list
  // }

  std::vector<int> targetMapProcessRanks;
  if (myRank == 0) {
    targetMapProcessRanks = {1, 2, 1, 2};
  }
  else if (myRank == 1) {
    targetMapProcessRanks = {2, 2};
  }
  // else if (myRank == 2) {
  // // do nothing; empty list
  // }

  const GST globalNumSourceMapGlobalIndices = 12;
  const GO indexBase = 0;
  return {comm, contiguous, localNumSourceMapGlobalIndices,
      sourceMapGlobalIndices,
      targetMapGlobalIndices,
      targetMapProcessRanks,
      globalNumSourceMapGlobalIndices,
      indexBase};
}

template<class LO, class GO, class NT>
void
runTest (Teuchos::FancyOStream& out,
         bool& success,
         const Issue2198TestInput<LO, GO, NT>& testInput,
         const std::string& testName)
{
  typedef ::Tpetra::Map<LO, GO, NT> map_type;
  typedef ::Tpetra::Import<LO, GO, NT> import_type;

  out << "Test " << testName << endl;
  Teuchos::OSTab tab1 (out);
  auto comm = testInput.comm;

  const bool lclTestInputLegit =
    testInput.targetMapGlobalIndices.size () ==
    testInput.targetMapProcessRanks.size ();
  const bool gblTestInputLegit = trueEverywhere (lclTestInputLegit, *comm);
  TEST_ASSERT( gblTestInputLegit );
  if (! gblTestInputLegit) {
    out << "Test input is broken; returning early" << std::endl;
    return;
  }

  auto sourceMap = makeSourceMapFromTestInput (testInput);
  auto expectedUnoptimizedTargetMap = makeUnoptimizedTargetMapFromTestInput (testInput, *sourceMap);
  RCP<import_type> expectedUnoptimizedImport =
    rcp (new import_type (sourceMap, expectedUnoptimizedTargetMap));

  RCP<const map_type> expectedOptimizedTargetMap;
  RCP<import_type> expectedOptimizedImport;
  {
    std::ostringstream errStrm;
    bool lclErr = false;
    using Tpetra::Details::makeOptimizedColMapAndImport;
    auto result = makeOptimizedColMapAndImport (out, lclErr, *sourceMap,
                                                *expectedUnoptimizedTargetMap,
                                                expectedUnoptimizedImport.getRawPtr (),
                                                true);
    const bool gblErr = ! trueEverywhere (! lclErr, *comm);
    TEST_ASSERT( ! gblErr );
    if (gblErr) {
      out << "makeOptimizedColMapAndImport FAILED; returning early!" << endl;
      return;
    }
    expectedOptimizedTargetMap = Teuchos::rcp (new map_type (result.first));
    expectedOptimizedImport = result.second;
  }

  out << "Source Map:" << endl;
  printMapCompactly (out, *sourceMap);

  out << "Expected unoptimized target Map:" << endl;
  printMapCompactly (out, *expectedUnoptimizedTargetMap);

  out << "Expected optimized target Map:" << endl;
  printMapCompactly (out, *expectedOptimizedTargetMap);

  for (bool optimized : {false, true}) {
    RCP<const import_type> expectedImport = optimized ?
      expectedOptimizedImport : expectedUnoptimizedImport;
    RCP<const map_type> expectedTargetMap = optimized ?
      expectedOptimizedTargetMap : expectedUnoptimizedTargetMap;

    std::unique_ptr<import_type> actualImport;
    std::ostringstream lclErrMsg;
    bool lclMadeItThrough = false;
    try {
      actualImport =
        std::unique_ptr<import_type> (new import_type (sourceMap,
                                                       testInput.targetMapGlobalIndices.data (),
                                                       testInput.targetMapProcessRanks.data (),
                                                       testInput.targetMapGlobalIndices.size (),
                                                       optimized));
      lclMadeItThrough = true;
    }
    catch (std::exception& e) {
      lclErrMsg << "Process " << comm->getRank ()
                << " threw an exception: " << e.what () << endl;
    }
    catch (...) {
      lclErrMsg << "Process " << comm->getRank ()
                << " threw an exception not a subclass of std::exception."
                << endl;
    }
    const bool gblMadeItThrough = trueEverywhere (lclMadeItThrough, *comm);
    TEST_ASSERT( gblMadeItThrough );
    if (! gblMadeItThrough) {
      Tpetra::Details::gathervPrint (out, lclErrMsg.str (), *comm);
    }
    else {
      // describe() is a collective, so we have to make sure that the
      // target Map is nonnull everywhere before calling describe().
      const bool gblTargetMapNonnull =
        trueEverywhere (! actualImport->getTargetMap ().is_null (), *comm);
      TEST_ASSERT( gblTargetMapNonnull );

      if (gblTargetMapNonnull) {
        out << "Actual unoptimized target Map:" << endl;
        printMapCompactly (out, * (actualImport->getTargetMap ()));

        const bool gblImportsSame =
          importsGloballySame (out, *expectedImport,
                               "expectedImport", *actualImport,
                               "actualImport");
        TEST_ASSERT( gblImportsSame );
      }
    }
  }
}


TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( ImportExport, Issue2198, LO, GO, NT )
{
  //typedef Tpetra::Map<LO, GO, NT> map_type; // unused
  //typedef Tpetra::Import<LO, GO, NT> import_type; // unused
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

  runTest (out, success, makeTest_A<LO, GO, NT> (comm), "A");
  runTest (out, success, makeTest_B<LO, GO, NT> (comm), "B");

  lclSuccess = success ? 1 : 0;
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY_CONST( gblSuccess, 1 );
}

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( ImportExport, Issue2198, LO, GO, NT )

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_LGN( UNIT_TEST_GROUP )
