// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Tpetra_ConfigDefs.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_Util.hpp> // sort2, merge2
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include "TpetraCore_ETIHelperMacros.h"

namespace { // anonymous

using Tpetra::global_size_t;
using Teuchos::Array;
using Teuchos::as;
using Teuchos::Comm;
using Teuchos::OSTab;
using Teuchos::ParameterList;
using Teuchos::ptr;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::REDUCE_MAX;
using Teuchos::REDUCE_MIN;
using Teuchos::reduceAll;
using std::endl;

const bool callFillComplete = true;
const bool tolerant = false;
// Whether to print copious debugging output to stderr when doing
// Matrix Market input and output.
const bool debug = false;

const char graph_symRealSmall[] =
"%%MatrixMarket matrix coordinate pattern general\n"
"5 5 13\n"
"1 1\n"
"1 2\n"
"2 1\n"
"2 2\n"
"2 3\n"
"3 2\n"
"3 3\n"
"3 4\n"
"4 3\n"
"4 4\n"
"4 5\n"
"5 4\n"
"5 5\n";


// Given an arbitrary Map, compute a Map containing all the GIDs
// in the same order as in (the one-to-one version of) map, but
// all owned exclusively by Proc 0.
template<class MapType>
Teuchos::RCP<const MapType>
computeGatherMap (Teuchos::RCP<const MapType> map,
                  const Teuchos::RCP<Teuchos::FancyOStream>& err=Teuchos::null,
                  const bool dbg=false)
{
  using Tpetra::createOneToOne;
  using Tpetra::global_size_t;
  using Teuchos::arcp;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::as;
  using Teuchos::Comm;
  using Teuchos::gather;
  using Teuchos::gatherv;
  using Teuchos::RCP;
  using std::endl;
  typedef typename MapType::local_ordinal_type LO;
  typedef typename MapType::global_ordinal_type GO;
  typedef typename MapType::node_type NT;

  RCP<const Comm<int> > comm = map->getComm ();
  const int numProcs = comm->getSize ();
  const int myRank = comm->getRank ();

  if (! err.is_null ()) {
    err->pushTab ();
  }
  if (dbg) {
    *err << myRank << ": computeGatherMap:" << endl;
  }
  if (! err.is_null ()) {
    err->pushTab ();
  }

  RCP<const MapType> oneToOneMap;
  if (map->isContiguous ()) {
    oneToOneMap = map; // contiguous Maps are always 1-to-1
  } else {
    if (dbg) {
      *err << myRank << ": computeGatherMap: Calling createOneToOne" << endl;
    }
    // It could be that Map is one-to-one, but the class doesn't
    // give us a way to test this, other than to create the
    // one-to-one Map.
    oneToOneMap = createOneToOne<LO, GO, NT> (map);
  }

  RCP<const MapType> gatherMap;
  if (numProcs == 1) {
    gatherMap = oneToOneMap;
  } else {
    if (dbg) {
      *err << myRank << ": computeGatherMap: Gathering Map counts" << endl;
    }
    // Gather each process' count of Map elements to Proc 0,
    // into the recvCounts array.  This will tell Proc 0 how
    // many GIDs to expect from each process when calling
    // MPI_Gatherv.  Counts and offsets are all int, because
    // that's what MPI uses.  Teuchos::as will at least prevent
    // bad casts to int in a dbg build.
    const int myEltCount = as<int> (oneToOneMap->getLocalNumElements ());
    Array<int> recvCounts (numProcs);
    const int rootProc = 0;
    gather<int, int> (&myEltCount, 1, recvCounts.getRawPtr (), 1, rootProc, *comm);

    ArrayView<const GO> myGlobalElts = oneToOneMap->getLocalElementList ();
    const int numMyGlobalElts = as<int> (myGlobalElts.size ());
    // Only Proc 0 needs to receive and store all the GIDs (from
    // all processes).
    ArrayRCP<GO> allGlobalElts;
    if (myRank == 0) {
      allGlobalElts = arcp<GO> (oneToOneMap->getGlobalNumElements ());
      std::fill (allGlobalElts.begin (), allGlobalElts.end (), 0);
    }

    if (dbg) {
      *err << myRank << ": computeGatherMap: Computing MPI_Gatherv "
        "displacements" << endl;
    }
    // Displacements for gatherv() in this case (gathering into a
    // contiguous array) are an exclusive partial sum (first entry is
    // zero, second starts the partial sum) of recvCounts.
    Array<int> displs (numProcs, 0);
    std::partial_sum (recvCounts.begin (), recvCounts.end () - 1,
                      displs.begin () + 1);
    if (dbg) {
      *err << myRank << ": computeGatherMap: Calling MPI_Gatherv" << endl;
    }
    gatherv<int, GO> (myGlobalElts.getRawPtr (), numMyGlobalElts,
                      allGlobalElts.getRawPtr (), recvCounts.getRawPtr (),
                      displs.getRawPtr (), rootProc, *comm);
    if (dbg) {
      *err << myRank << ": computeGatherMap: Creating gather Map" << endl;
    }
    // Create a Map with all the GIDs, in the same order as in
    // the one-to-one Map, but owned by Proc 0.
    ArrayView<const GO> allElts (NULL, 0);
    if (myRank == 0) {
      allElts = allGlobalElts ();
    }
    const global_size_t INVALID = Teuchos::OrdinalTraits<global_size_t>::invalid ();
    gatherMap = rcp (new MapType (INVALID, allElts,
                                  oneToOneMap->getIndexBase (),
                                  comm));
  }
  if (! err.is_null ()) {
    err->popTab ();
  }
  if (dbg) {
    *err << myRank << ": computeGatherMap: done" << endl;
  }
  if (! err.is_null ()) {
    err->popTab ();
  }
  return gatherMap;
}

// Input matrices must be fill complete.
template<class CrsGraphType>
bool
compareCrsGraphMaps (const CrsGraphType& A_orig, const CrsGraphType& A, Teuchos::FancyOStream& out)
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::Comm;
  using Teuchos::RCP;
  using Teuchos::reduceAll;
  using Teuchos::REDUCE_MIN;

  Teuchos::OSTab tab (Teuchos::rcpFromRef (out));

  bool globalAllSame = true;
  if (! A_orig.getRowMap ()->isSameAs (* (A.getRowMap ()))) {
    out << "Row Maps are not the same" << endl;
    globalAllSame = false;
  }
  if (! A_orig.getColMap ()->isSameAs (* (A.getColMap ()))) {
    out << "Column Maps are not the same" << endl;
    globalAllSame = false;
  }
  if (! A_orig.getDomainMap ()->isSameAs (* (A.getDomainMap ()))) {
    out << "Domain Maps are not the same" << endl;
    globalAllSame = false;
  }
  if (! A_orig.getRangeMap ()->isSameAs (* (A.getRangeMap ()))) {
    out << "Range Maps are not the same" << endl;
    globalAllSame = false;
  }
  if (globalAllSame) {
    out << "All Maps are the same" << endl;
  }
  return globalAllSame;
}

// Input matrices must be fill complete, and all four of their Maps
// (row, column, domain, and range) must be the same.
template<class CrsGraphType>
bool
compareCrsGraph (const CrsGraphType& A_orig, const CrsGraphType& A, Teuchos::FancyOStream& out)
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::Comm;
  using Teuchos::RCP;
  using Teuchos::reduceAll;
  using Teuchos::REDUCE_MIN;
  typedef typename CrsGraphType::global_ordinal_type GO;
  typedef typename ArrayView<const GO>::size_type size_type;
  typedef typename CrsGraphType::nonconst_global_inds_host_view_type gids_type;

  Teuchos::OSTab tab (Teuchos::rcpFromRef (out));
  int localEqual = 1;

  //
  // Are my local matrices equal?
  //
  gids_type indOrig, ind;
  size_t numEntriesOrig = 0;
  size_t numEntries = 0;

  ArrayView<const GO> localElts = A.getRowMap ()->getLocalElementList ();
  const size_type numLocalElts = localElts.size ();
  for (size_type i = 0; i < numLocalElts; ++i) {
    const GO globalRow = localElts[i];
    numEntriesOrig = A_orig.getNumEntriesInGlobalRow (globalRow);
    numEntries = A.getNumEntriesInGlobalRow (globalRow);

    if (numEntriesOrig != numEntries) {
      localEqual = 0;
      break;
    }
    Kokkos::resize(indOrig,numEntriesOrig);
    A_orig.getGlobalRowCopy (globalRow, indOrig, numEntriesOrig);
    Kokkos::resize(ind,numEntries);
    A.getGlobalRowCopy (globalRow, ind, numEntries);

    // Global row entries are not necessarily sorted.  Sort them so
    // we can compare them.
    Tpetra::sort (indOrig, indOrig.extent(0));
    Tpetra::sort (ind, ind.extent(0));

    for (size_t k = 0; k < numEntries; ++k) {
      // Indices should be _exactly_ equal.
      if (indOrig[k] != ind[k]) {
        localEqual = 0;
        break;
      }
    }
  }

  RCP<const Comm<int> > comm = A.getRowMap ()->getComm ();
  int globalEqual = 0;
  reduceAll<int, int> (*comm, REDUCE_MIN, 1, &localEqual, &globalEqual);
  return globalEqual == 1;
}


template<class LocalOrdinalType, class GlobalOrdinalType, class NodeType>
Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinalType, GlobalOrdinalType, NodeType> >
createSymRealSmall (const Teuchos::RCP<const Tpetra::Map<LocalOrdinalType, GlobalOrdinalType, NodeType> >& rowMap,
                    Teuchos::FancyOStream& out,
                    const bool dbg)
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  typedef LocalOrdinalType LO;
  typedef GlobalOrdinalType GO;
  typedef NodeType NT;
  typedef Tpetra::global_size_t GST;
  typedef Tpetra::Map<LO, GO, NT> map_type;
  typedef Tpetra::Export<LO, GO, NT> export_type;
  typedef Tpetra::CrsGraph<LO, GO, NT> graph_type;

  RCP<const Teuchos::Comm<int> > comm = rowMap->getComm ();
  const int myRank = comm->getRank ();

  const GST globalNumElts = rowMap->getGlobalNumElements ();
  const size_t myNumElts = (myRank == 0) ? as<size_t> (globalNumElts) : 0;
  RCP<const map_type> gatherRowMap = computeGatherMap (rowMap, rcpFromRef (out), dbg);
  graph_type A_gather (gatherRowMap, 3);

  if (myRank == 0) {
    Array<GO> ind (3);
    for (size_t myRow = 0; myRow < myNumElts; ++myRow) {
      const GO globalRow = gatherRowMap->getGlobalElement (myRow);
      if (globalRow == gatherRowMap->getMinAllGlobalIndex ()) {
        ind[0] = globalRow;
        ind[1] = globalRow + 1;
        A_gather.insertGlobalIndices (globalRow, ind.view (0, 2));
      }
      else if (globalRow == gatherRowMap->getMaxAllGlobalIndex ()) {
        ind[0] = globalRow - 1;
        ind[1] = globalRow;
        A_gather.insertGlobalIndices (globalRow, ind.view (0, 2));
      }
      else {
        ind[0] = globalRow - 1;
        ind[1] = globalRow;
        ind[2] = globalRow + 1;
        A_gather.insertGlobalIndices (globalRow, ind.view (0, 3));
      }
    }
  }
  A_gather.fillComplete (rowMap, rowMap);
  RCP<graph_type> A = rcp (new graph_type (rowMap, as<size_t> (0)));
  export_type exp (gatherRowMap, rowMap);
  A->doExport (A_gather, exp, Tpetra::INSERT);
  A->fillComplete (rowMap, rowMap);
  return A;
}

template<class LocalOrdinalType, class GlobalOrdinalType, class NodeType>
bool
testCrsGraph (Teuchos::FancyOStream& out, const GlobalOrdinalType indexBase)
{
  typedef LocalOrdinalType LO;
  typedef GlobalOrdinalType GO;
  typedef NodeType NT;
  typedef Tpetra::Map<LO, GO, NT> map_type;
  typedef Tpetra::CrsGraph<LO, GO, NT> crs_graph_type;
  typedef Tpetra::CrsMatrix<>::scalar_type scalar_type;
  typedef Tpetra::CrsMatrix<scalar_type,LO, GO, NT> crs_matrix_type;
  bool result = true; // current Boolean result; reused below
  bool success = true;

  out << "Test: CrsGraph Matrix Market I/O, w/ Map with index base "
      << indexBase << endl;
  OSTab tab1 (out);

  RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();

  out << "Original sparse graph:" << endl;
  out << graph_symRealSmall << endl;

  out << "Creating the row Map" << endl;
  const global_size_t globalNumElts = 5;
  RCP<const map_type> rowMap =
    rcp (new map_type (globalNumElts, indexBase, comm,
                       Tpetra::GloballyDistributed));

  out << "Reading in the graph" << endl;
  std::istringstream inStr (graph_symRealSmall);
  RCP<const map_type> colMap;
  RCP<const map_type> domainMap = rowMap;
  RCP<const map_type> rangeMap = rowMap;
  typedef Tpetra::MatrixMarket::Reader<crs_matrix_type> reader_type;
  out << "Instantiating the graph" << endl;
  RCP<crs_graph_type> A =
    reader_type::readSparseGraph (inStr, rowMap, colMap, domainMap, rangeMap,
                             callFillComplete, tolerant, debug);

  out << "Creating original graph" << endl;
  RCP<crs_graph_type> A_orig =
    createSymRealSmall<LO, GO, NT> (rowMap, out, debug);

  out << "Comparing read-in graph to original graph" << endl;
  result = compareCrsGraph<crs_graph_type> (*A_orig, *A, out);
  bool local_success = true;
  TEUCHOS_TEST_EQUALITY( result, true, out, local_success );
  success = success && local_success;

  out << "Writing out the original graph" << endl;
  std::ostringstream outStr;
  typedef Tpetra::MatrixMarket::Writer<crs_matrix_type> writer_type;
  writer_type::writeSparseGraph (outStr, *A_orig, debug);

  out << "Reading it in again and comparing with original graph" << endl;
  std::istringstream inStr2 (outStr.str ());
  RCP<crs_graph_type> A_orig2 =
    reader_type::readSparseGraph (inStr2, rowMap, colMap, domainMap, rangeMap,
                             callFillComplete, tolerant, debug);
  result = compareCrsGraph<crs_graph_type> (*A_orig, *A_orig2, out);
  local_success = true;
  TEUCHOS_TEST_EQUALITY( result, true, out, local_success );
  success = success && local_success;

  // mfh 01 Jul 2018: We just need this to compile for now.  Later, we
  // should go back and add a real test.  Note that this interface
  // doesn't necessarily promise the same Maps as the interface above.
  out << "Test fix for #3043" << endl;
  {
    Teuchos::RCP<Teuchos::ParameterList> fakeParams;
    {
      std::istringstream istr3 (outStr.str ());
      RCP<crs_graph_type> A_orig3 =
        reader_type::readSparseGraph (istr3,
                                      rowMap->getComm (),
                                      callFillComplete, tolerant, debug);
      // result = compareCrsGraph<crs_graph_type> (*A_orig, *A_orig3, out);
      // local_success = true;
      local_success = A_orig3.get () != nullptr;
      TEUCHOS_TEST_EQUALITY( result, true, out, local_success );
      success = success && local_success;
    }
    {
      std::istringstream istr4 (outStr.str ());
      RCP<crs_graph_type> A_orig4 =
        reader_type::readSparseGraph (istr4,
                                      rowMap->getComm (),
                                      fakeParams,
                                      fakeParams,
                                      tolerant, debug);
      // result = compareCrsGraph<crs_graph_type> (*A_orig, *A_orig4, out);
      // local_success = true;
      local_success = A_orig4.get () != nullptr;
      TEUCHOS_TEST_EQUALITY( result, true, out, local_success );
      success = success && local_success;
    }
  }
  return success;
}

} // namespace (anonymous)

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraphOutputInput, IndexBase0, LO, GO, NT )
{
  const GO indexBase = 0;
  success = testCrsGraph<LO, GO, NT> (out, indexBase);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( CrsGraphOutputInput, IndexBase1, LO, GO, NT )
{
  const GO indexBase = 1;
  success = testCrsGraph<LO, GO, NT> (out, indexBase);
}

// We instantiate tests for all combinations of the following parameters:
//   - indexBase = {0, 1}
//   - Scalar = {double, float}

#  define UNIT_TEST_GROUP( LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraphOutputInput, IndexBase0, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( CrsGraphOutputInput, IndexBase1, LO, GO, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_LGN( UNIT_TEST_GROUP )

