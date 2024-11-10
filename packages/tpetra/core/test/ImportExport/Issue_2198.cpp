// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Test for Github Issue #2198 (feature request for new Import constructor).
//
// mfh 10 May 2018: See GitHub Issue #2564 for why we split this test
// into two separate cases: 3 MPI processes, and 5 MPI processes.

namespace { // (anonymous)
constexpr int minNumProcsForTestA = 5;
constexpr int minNumProcsForTestB = 3;
} // namespace (anonymous)

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
using GST = Tpetra::global_size_t;

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
  if (map.getLocalNumElements () != 0) {
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

  {
    auto X_permuteFromLIDs = X.getPermuteFromLIDs_dv ();
    auto Y_permuteFromLIDs = Y.getPermuteFromLIDs_dv ();
    if (X_permuteFromLIDs.extent (0) != Y_permuteFromLIDs.extent (0)) {
      out << prefix << X_name << ".getPermuteFromLIDs_dv().extent(0)="
	  << X_permuteFromLIDs.extent (0)
	  << " != " << Y_name << ".getPermuteFromLIDs_dv().extent(0)="
	  << Y_permuteFromLIDs.extent (0)
	  << endl;
      same = false;
    }
    else {
      auto X_ptr = X_permuteFromLIDs.view_host ().data ();
      const auto size = X_permuteFromLIDs.view_host ().extent (0);
      auto Y_ptr = Y_permuteFromLIDs.view_host ().data ();
      
      if (! std::equal (X_ptr, X_ptr + size, Y_ptr)) {
	out << prefix << X_name << ".getPermuteFromLIDs_dv().view_host()"
	    << " != " << Y_name << ".getPermuteFromLIDs_dv().view_host()"
	    << endl;
	same = false;
      }
    }
  }

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

  {
    auto X_permuteToLIDs = X.getPermuteToLIDs_dv ();
    auto Y_permuteToLIDs = Y.getPermuteToLIDs_dv ();
    if (X_permuteToLIDs.extent (0) != Y_permuteToLIDs.extent (0)) {
      out << prefix << X_name << ".getPermuteToLIDs_dv().extent(0)="
	  << X_permuteToLIDs.extent (0)
	  << " != " << Y_name << ".getPermuteToLIDs_dv().extent(0)="
	  << Y_permuteToLIDs.extent (0)
	  << endl;
      same = false;
    }
    else {
      auto X_ptr = X_permuteToLIDs.view_host ().data ();
      const auto size = X_permuteToLIDs.view_host ().extent (0);
      auto Y_ptr = Y_permuteToLIDs.view_host ().data ();
      
      if (! std::equal (X_ptr, X_ptr + size, Y_ptr)) {
	out << prefix << X_name << ".getPermuteToLIDs_dv().view_host()"
	    << " != " << Y_name << ".getPermuteToLIDs_dv().view_host()"
	    << endl;
	same = false;
      }
    }
  }

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

  {
    auto X_exportLIDs = X.getExportLIDs_dv ();
    auto Y_exportLIDs = Y.getExportLIDs_dv ();
    if (X_exportLIDs.extent (0) != Y_exportLIDs.extent (0)) {
      out << prefix << X_name << ".getExportLIDs_dv().extent(0)="
	  << X_exportLIDs.extent (0)
	  << " != " << Y_name << ".getExportLIDs_dv().extent(0)="
	  << Y_exportLIDs.extent (0)
	  << endl;
      same = false;
    }
    else {
      auto X_ptr = X_exportLIDs.view_host ().data ();
      const auto size = X_exportLIDs.view_host ().extent (0);
      auto Y_ptr = Y_exportLIDs.view_host ().data ();
      
      if (! std::equal (X_ptr, X_ptr + size, Y_ptr)) {
	out << prefix << X_name << ".getExportLIDs_dv().view_host()"
	    << " != " << Y_name << ".getExportLIDs_dv().view_host()"
	    << endl;
	same = false;
      }
    }
  }

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

  {
    auto X_remoteLIDs = X.getRemoteLIDs_dv ();
    auto Y_remoteLIDs = Y.getRemoteLIDs_dv ();
    if (X_remoteLIDs.extent (0) != Y_remoteLIDs.extent (0)) {
      out << prefix << X_name << ".getRemoteLIDs_dv().extent(0)="
	  << X_remoteLIDs.extent (0)
	  << " != " << Y_name << ".getRemoteLIDs_dv().extent(0)="
	  << Y_remoteLIDs.extent (0)
	  << endl;
      same = false;
    }
    else {
      auto X_ptr = X_remoteLIDs.view_host ().data ();
      const auto size = X_remoteLIDs.view_host ().extent (0);
      auto Y_ptr = Y_remoteLIDs.view_host ().data ();
      
      if (! std::equal (X_ptr, X_ptr + size, Y_ptr)) {
	out << prefix << X_name << ".getRemoteLIDs_dv().view_host()"
	    << " != " << Y_name << ".getRemoteLIDs_dv().view_host()"
	    << endl;
	same = false;
      }
    }
  }

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
  std::vector<GO> remoteGlobalIndices;
  std::vector<GO> optimizedRemoteGlobalIndices;
  std::vector<int> remoteProcessRanks;
  std::vector<int> optimizedRemoteProcessRanks;
  const GST globalNumSourceMapGlobalIndices;
  const GO indexBase;
};

template<class LO, class GO, class NT>
Teuchos::RCP<const Tpetra::Map<LO, GO, NT> >
makeSourceMapFromTestInput (const Issue2198TestInput<LO, GO, NT>& test)
{
  using map_type = Tpetra::Map<LO, GO, NT>;
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
makeTargetMapFromTestInput (const Issue2198TestInput<LO, GO, NT>& testInput,
                            const ::Tpetra::Map<LO, GO, NT>& sourceMap,
                            const bool optimized)
{
  typedef Tpetra::Map<LO, GO, NT> map_type;

  const std::vector<GO>& remoteGlobalIndices = optimized ?
    testInput.optimizedRemoteGlobalIndices : testInput.remoteGlobalIndices;

  const LO numLclSrcGids = static_cast<LO> (sourceMap.getLocalNumElements ());
  const LO numInputGids = static_cast<LO> (remoteGlobalIndices.size ());
  const LO numLclTgtGids = numLclSrcGids + numInputGids;
  std::vector<GO> tgtGids (numLclTgtGids);

  if (sourceMap.isContiguous ()) {
    GO curTgtGid = sourceMap.getMinGlobalIndex ();
    for (LO k = 0; k < numLclSrcGids; ++k, ++curTgtGid) {
      tgtGids[k] = curTgtGid;
    }
  }
  else {
    auto srcGids = sourceMap.getLocalElementList ();
    for (LO k = 0; k < numLclSrcGids; ++k) {
      tgtGids[k] = srcGids[k];
    }
  }
  std::copy (remoteGlobalIndices.begin (),
             remoteGlobalIndices.end (),
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
  const int numProcs = comm->getSize ();
  TEUCHOS_TEST_FOR_EXCEPTION
    (numProcs < minNumProcsForTestA, std::invalid_argument, "Test A needs at "
     "least " << minNumProcsForTestA << "MPI processes, but the given "
     "communicator has only " << numProcs
     << " process" << (numProcs != 1 ? "es" : "") << ".");

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

  std::vector<GO> remoteGlobalIndices;
  std::vector<GO> optimizedRemoteGlobalIndices;
  if (myRank == 0) {
    remoteGlobalIndices = {1, 3, 7};
    optimizedRemoteGlobalIndices = {1, 3, 7};
  }
  else if (myRank == 1) {
    remoteGlobalIndices = {3, 9, 10};
    optimizedRemoteGlobalIndices = {10, 3, 9};
  }
  else if (myRank == 2) {
    remoteGlobalIndices = {6, 7, 13};
    optimizedRemoteGlobalIndices = {6, 7, 13};
  }
  else if (myRank == 3) {
    remoteGlobalIndices = {9, 13, 15};
    optimizedRemoteGlobalIndices = {15, 9, 13};
  }
  else if (myRank == 4) {
    remoteGlobalIndices = {4, 5, 11, 12};
    optimizedRemoteGlobalIndices = {4, 5, 11, 12};
  }

  std::vector<int> remoteProcessRanks;
  std::vector<int> optimizedRemoteProcessRanks;
  if (myRank == 0) {
    remoteProcessRanks = {1, 4, 4};
    optimizedRemoteProcessRanks = {1, 4, 4};
  }
  else if (myRank == 1) {
    remoteProcessRanks = {4, 4, 3};
    optimizedRemoteProcessRanks = {3, 4, 4};
  }
  else if (myRank == 2) {
    remoteProcessRanks = {0, 4, 4};
    optimizedRemoteProcessRanks = {0, 4, 4};
  }
  else if (myRank == 3) {
    remoteProcessRanks = {4, 4, 2};
    optimizedRemoteProcessRanks = {2, 4, 4};
  }
  else if (myRank == 4) {
    remoteProcessRanks = {0, 1, 2, 3};
    optimizedRemoteProcessRanks = {0, 1, 2, 3};
  }

  const GST globalNumSourceMapGlobalIndices = 17;
  const GO indexBase = 0;
  return {comm,
      contiguous,
      localNumSourceMapGlobalIndices,
      sourceMapGlobalIndices,
      remoteGlobalIndices,
      optimizedRemoteGlobalIndices,
      remoteProcessRanks,
      optimizedRemoteProcessRanks,
      globalNumSourceMapGlobalIndices,
      indexBase};
}

template<class LO, class GO, class NT>
Issue2198TestInput<LO, GO, NT>
makeTest_B (const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
{
  const int numProcs = comm->getSize ();
  TEUCHOS_TEST_FOR_EXCEPTION
    (numProcs < minNumProcsForTestB, std::invalid_argument, "Test A needs at "
     "least " << minNumProcsForTestB << "MPI processes, but the given "
     "communicator has only " << numProcs
     << " process" << (numProcs != 1 ? "es" : "") << ".");

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

  std::vector<GO> remoteGlobalIndices;
  std::vector<GO> optimizedRemoteGlobalIndices;
  if (myRank == 0) {
    remoteGlobalIndices = {5, 11, 3, 8};
    optimizedRemoteGlobalIndices = {5, 3, 11, 8};
  }
  else if (myRank == 1) {
    remoteGlobalIndices = {10, 7};
    optimizedRemoteGlobalIndices = {10, 7};
  }
  // else if (myRank == 2) {
  // // do nothing; empty list
  // }

  std::vector<int> remoteProcessRanks;
  std::vector<int> optimizedRemoteProcessRanks;
  if (myRank == 0) {
    remoteProcessRanks = {1, 2, 1, 2};
    optimizedRemoteProcessRanks = {1, 1, 2, 2};
  }
  else if (myRank == 1) {
    remoteProcessRanks = {2, 2};
    optimizedRemoteProcessRanks = {2, 2};
  }
  // else if (myRank == 2) {
  // // do nothing; empty list
  // }

  const GST globalNumSourceMapGlobalIndices = 12;
  const GO indexBase = 0;
  return {comm, contiguous, localNumSourceMapGlobalIndices,
      sourceMapGlobalIndices,
      remoteGlobalIndices,
      optimizedRemoteGlobalIndices,
      remoteProcessRanks,
      optimizedRemoteProcessRanks,
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
  using map_type = ::Tpetra::Map<LO, GO, NT>;
  using import_type = ::Tpetra::Import<LO, GO, NT>;

  out << "Test " << testName << endl;
  Teuchos::OSTab tab1 (out);
  auto comm = testInput.comm;

  const bool lclTestInputLegit =
    (testInput.remoteGlobalIndices.size () ==
     testInput.remoteProcessRanks.size ()) &&
    (testInput.optimizedRemoteGlobalIndices.size () ==
     testInput.optimizedRemoteProcessRanks.size ());
  const bool gblTestInputLegit = trueEverywhere (lclTestInputLegit, *comm);
  TEST_ASSERT( gblTestInputLegit );
  if (! gblTestInputLegit) {
    out << "Test input is broken; returning early" << std::endl;
    return;
  }

  auto sourceMap = makeSourceMapFromTestInput (testInput);
  auto expUnoptTgtMap = makeTargetMapFromTestInput (testInput, *sourceMap, false);
  RCP<import_type> expUnoptImport (new import_type (sourceMap, expUnoptTgtMap));

  auto expOptTgtMap = makeTargetMapFromTestInput (testInput, *sourceMap, true);
  RCP<import_type> expOptImport (new import_type (sourceMap, expOptTgtMap));

  out << "Source Map:" << endl;
  printMapCompactly (out, *sourceMap);
  out << endl;

  for (bool optimized : {false, true}) {
    out << "Optimized: " << (optimized ? "true" : "false") << endl;
    Teuchos::OSTab tab2 (out);

    RCP<const import_type> expImport = optimized ? expOptImport : expUnoptImport;
    RCP<const map_type> expTgtMap = optimized ? expOptTgtMap : expUnoptTgtMap;

    out << "Expected target Map:" << endl;
    printMapCompactly (out, *expTgtMap);

    std::unique_ptr<import_type> actualImport;
    std::ostringstream lclErrMsg;
    bool lclMadeItThrough = false;
    try {
      actualImport =
        std::unique_ptr<import_type> (new import_type (sourceMap,
                                                       testInput.remoteGlobalIndices.data (),
                                                       testInput.remoteProcessRanks.data (),
                                                       testInput.remoteGlobalIndices.size (),
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
        out << "Actual target Map:" << endl;
        printMapCompactly (out, * (actualImport->getTargetMap ()));

        const bool gblImportsSame =
          importsGloballySame (out, *expImport,
                               "expImport", *actualImport,
                               "actualImport");
        TEST_ASSERT( gblImportsSame );
      }
    }

    const bool allSuccessfulThusFar = trueEverywhere (success, *comm);
    if (! allSuccessfulThusFar) {
      return;
    }
  } // for optimized in {false, true}


  // Now test against results of makeOptimizedColMap.
  {
    out << "Test against makeOptimizedColMap" << endl;
    Teuchos::OSTab tab2 (out);

    std::ostringstream errStrm;
    bool lclErr = false;
    Teuchos::RCP<const map_type> actualTgtMap =
      Tpetra::Details::makeOptimizedColMap (out, lclErr, *sourceMap,
                                            *expUnoptTgtMap,
                                            expUnoptImport.getRawPtr ());
    const bool gblMadeItThrough = trueEverywhere (! lclErr, *comm);
    TEST_ASSERT( gblMadeItThrough );
    if (! gblMadeItThrough) {
      out << "makeOptimizedColMap FAILED; returning early!" << endl;
      return;
    }

    out << "Actual target Map:" << endl;
    printMapCompactly (out, *actualTgtMap);

    const import_type actualImport (sourceMap, actualTgtMap);
    const bool gblImportsSame =
      importsGloballySame (out, *expOptImport, "expOptImport",
                           actualImport, "actualImport");
    TEST_ASSERT( gblImportsSame );
  }
}


TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( ImportExport, Issue2198, LO, GO, NT )
{
  int lclSuccess = 1;
  int gblSuccess = 1;

  out << "Test for Issue 2198 (new Import constructor)" << endl;
  Teuchos::OSTab tab0 (out);

  auto comm = Teuchos::DefaultComm<int>::getComm ();
  const int numProcs = comm->getSize ();
  constexpr int minNumProcs = minNumProcsForTestA < minNumProcsForTestB ?
    minNumProcsForTestA : minNumProcsForTestB;
  if (numProcs < minNumProcs) {
    out << "Test FAILED; test requires at least " << minNumProcs << " MPI "
        << "processes, but the input communicator has " << numProcs
        << "process" << (numProcs != 1 ? "es" : "") << "." << endl;
    success = false;
    return;
  }

  bool didSomeTestRun = false;

  // mfh 10 May 2018: See GitHub Issue #2564 for why we split this
  // test into two separate cases, depending on the number of MPI
  // processes that the user gives us.
  if (numProcs >= minNumProcsForTestA && ! didSomeTestRun) {
    runTest (out, success, makeTest_A<LO, GO, NT> (comm), "A");
    didSomeTestRun = true;
  }
  if (numProcs >= minNumProcsForTestB && ! didSomeTestRun) {
    runTest (out, success, makeTest_B<LO, GO, NT> (comm), "B");
    didSomeTestRun = true;
  }

  // Don't let the test get by with doing nothing.
  if (! didSomeTestRun) {
    out << "Test FAILED, since no tests ran!" << endl;
    success = false;
    return;
  }

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

using default_local_ordinal_type = Tpetra::Map<>::local_ordinal_type;
using default_global_ordinal_type = Tpetra::Map<>::global_ordinal_type;
using default_node_type = Tpetra::Map<>::node_type;

UNIT_TEST_GROUP( default_local_ordinal_type, default_global_ordinal_type, default_node_type )

//TPETRA_INSTANTIATE_LGN( UNIT_TEST_GROUP )
