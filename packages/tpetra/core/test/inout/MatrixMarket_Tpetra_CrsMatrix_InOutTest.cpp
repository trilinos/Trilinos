// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_Util.hpp> // sort2, merge2
#include <Tpetra_TestingUtilities.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include "TpetraCore_ETIHelperMacros.h"

namespace { // anonymous

using Tpetra::global_size_t;
using Tpetra::TestingUtilities::arcp_from_view;
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

const char matrix_symRealSmall[] =
"%%MatrixMarket matrix coordinate real general\n"
"5 5 13\n"
"1 1  2.0000000000000e+00\n"
"1 2  -1.0000000000000e+00\n"
"2 1  -1.0000000000000e+00\n"
"2 2  2.0000000000000e+00\n"
"2 3  -1.0000000000000e+00\n"
"3 2  -1.0000000000000e+00\n"
"3 3  2.0000000000000e+00\n"
"3 4  -1.0000000000000e+00\n"
"4 3  -1.0000000000000e+00\n"
"4 4  2.0000000000000e+00\n"
"4 5  -1.0000000000000e+00\n"
"5 4  -1.0000000000000e+00\n"
"5 5  2.0000000000000e+00\n";


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
template<class CrsMatrixType>
bool
compareCrsMatrixMaps (const CrsMatrixType& A_orig, const CrsMatrixType& A, Teuchos::FancyOStream& out)
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::Comm;
  using Teuchos::RCP;
  using Teuchos::reduceAll;
  using Teuchos::REDUCE_MIN;
  // typedef typename CrsMatrixType::scalar_type ST;
  // typedef typename CrsMatrixType::global_ordinal_type GO;
  // typedef typename ArrayView<const GO>::size_type size_type;

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
template<class CrsMatrixType>
bool
compareCrsMatrix (const CrsMatrixType& A_orig, const CrsMatrixType& A, Teuchos::FancyOStream& out)
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::Comm;
  using Teuchos::RCP;
  using Teuchos::reduceAll;
  using Teuchos::REDUCE_MIN;
  typedef typename CrsMatrixType::global_ordinal_type GO;
  typedef typename ArrayView<const GO>::size_type size_type;
  typedef typename CrsMatrixType::nonconst_global_inds_host_view_type gids_type;
  typedef typename CrsMatrixType::nonconst_values_host_view_type vals_type;

  Teuchos::OSTab tab (Teuchos::rcpFromRef (out));
  int localEqual = 1;

  //
  // Are my local matrices equal?
  //
  gids_type indOrig, ind;
  vals_type valOrig, val;
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
    Kokkos::resize(valOrig,numEntriesOrig);
    A_orig.getGlobalRowCopy (globalRow, indOrig, valOrig, numEntriesOrig);
    Kokkos::resize(ind,numEntries);
    Kokkos::resize(val,numEntries);
    A.getGlobalRowCopy (globalRow, ind, val, numEntries);

    // Global row entries are not necessarily sorted.  Sort them so
    Tpetra::sort2 (indOrig, indOrig.extent(0), valOrig);
    Tpetra::sort2 (ind, ind.extent(0), val);

    for (size_t k = 0; k < numEntries; ++k) {
      // Values should be _exactly_ equal.
      if (indOrig[k] != ind[k] || valOrig[k] != val[k]) {
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


// Input matrices must be fill complete, and all four of their Maps
// (row, column, domain, and range) must be the same.
template<class CrsMatrixType>
bool
compareCrsMatrixValues (const CrsMatrixType& A_orig,
                        const CrsMatrixType& A,
                        Teuchos::FancyOStream& out)
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::ArrayRCP;
  using Teuchos::Comm;
  using Teuchos::RCP;
  using Teuchos::reduceAll;
  using Teuchos::REDUCE_MIN;
  using std::endl;
  typedef typename CrsMatrixType::scalar_type ST;
  typedef typename CrsMatrixType::global_ordinal_type GO;
  typedef typename ArrayView<const GO>::size_type size_type;
  typedef Teuchos::ScalarTraits<ST> STS;
  typedef typename STS::magnitudeType MT;
  typedef Teuchos::ScalarTraits<MT> STM;
  typedef typename CrsMatrixType::nonconst_global_inds_host_view_type gids_type;
  typedef typename CrsMatrixType::nonconst_values_host_view_type vals_type;

  Teuchos::OSTab tab (Teuchos::rcpFromRef (out));

  //
  // Are my local matrices equal?
  //
 //
  gids_type indOrig_v, ind_v;
  vals_type valOrig_v, val_v;
  size_t numEntriesOrig = 0;
  size_t numEntries = 0;

  ArrayView<const GO> localElts = A.getRowMap ()->getLocalElementList ();
  const size_type numLocalElts = localElts.size ();
  MT localDiff = STM::zero (); // \sum_{i,j} |A_orig(i,j) - A(i,j)| locally
  for (size_type i = 0; i < numLocalElts; ++i) {
    const GO globalRow = localElts[i];
    numEntriesOrig = A_orig.getNumEntriesInGlobalRow (globalRow);
    numEntries = A.getNumEntriesInGlobalRow (globalRow);

    Kokkos::resize(indOrig_v,numEntriesOrig);
    Kokkos::resize(valOrig_v,numEntriesOrig);
    A_orig.getGlobalRowCopy (globalRow, indOrig_v, valOrig_v, numEntriesOrig);
    Kokkos::resize(ind_v,numEntries);
    Kokkos::resize(val_v,numEntries);
    A.getGlobalRowCopy (globalRow, ind_v , val_v, numEntries);

    // Global row entries are not necessarily sorted.  Sort them
    // (and their values with them) so we can merge their values.
    Tpetra::sort2 (indOrig_v, indOrig_v.extent(0), valOrig_v);
    Tpetra::sort2 (ind_v, ind_v.extent(0), val_v);

    auto indOrig = arcp_from_view(indOrig_v);
    auto ind = arcp_from_view(ind_v);
    auto valOrig = arcp_from_view(valOrig_v);
    auto val = arcp_from_view(val_v);

    //
    // Merge repeated values in each set of indices and values.
    //

    typename ArrayRCP<GO>::iterator indOrigIter = indOrig.begin ();
    typename ArrayRCP<ST>::iterator valOrigIter = valOrig.begin ();
    typename ArrayRCP<GO>::iterator indOrigEnd = indOrig.end ();
    typename ArrayRCP<ST>::iterator valOrigEnd = valOrig.end ();
    Tpetra::merge2 (indOrigEnd, valOrigEnd, indOrigIter, indOrigEnd, valOrigIter, valOrigEnd);

  
    typename ArrayRCP<GO>::iterator indIter = ind.begin ();
    typename ArrayRCP<ST>::iterator valIter = val.begin ();
    typename ArrayRCP<GO>::iterator indEnd = ind.end ();
    typename ArrayRCP<ST>::iterator valEnd = val.end ();
    Tpetra::merge2 (indEnd, valEnd, indIter, indEnd, valIter, valEnd);

    //
    // Compare the merged sets of entries.
    //

    indOrigIter = indOrig.begin ();
    indIter = ind.begin ();
    valOrigIter = valOrig.begin ();
    valIter = val.begin ();
    while (indOrigIter != indOrigEnd && indIter != indEnd) {
      const GO j_orig = *indOrigIter;
      const GO j = *indIter;

      if (j_orig < j) { // entry is in A_orig, not in A
        localDiff += STS::magnitude (*valOrigIter);
        ++indOrigIter;
        ++valOrigIter;
      } else if (j_orig > j) { // entry is in A, not A_orig
        localDiff += STS::magnitude (*valIter);
        ++indIter;
        ++valIter;
      } else { // j_orig == j: entry is in both matrices
        localDiff += STS::magnitude (*valOrigIter - *valIter);
        ++indOrigIter;
        ++valOrigIter;
        ++indIter;
        ++valIter;
      }
    }
  }

  RCP<const Comm<int> > comm = A.getRowMap ()->getComm ();
  const int myRank = comm->getRank ();
  std::ostringstream os;
  os << "Values are ";
  if (localDiff == STM::zero ()) {
    os << "the same on process " << myRank << endl;
  } else {
    os << "NOT the same on process " << myRank
       << ": \\sum_{i,j} |A_orig(i,j) - A(i,j)| = " << localDiff << endl;
  }
  out << os.str ();

  int globalEqual = 0;
  int localEqual = (localDiff == STM::zero ()) ? 1 : 0;
  reduceAll<int, int> (*comm, REDUCE_MIN, 1, &localEqual, &globalEqual);
  return globalEqual == 1;
}

template<class ScalarType, class LocalOrdinalType, class GlobalOrdinalType, class NodeType>
Teuchos::RCP<Tpetra::CrsMatrix<ScalarType, LocalOrdinalType, GlobalOrdinalType, NodeType> >
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
  typedef ScalarType ST;
  typedef LocalOrdinalType LO;
  typedef GlobalOrdinalType GO;
  typedef NodeType NT;
  typedef Tpetra::global_size_t GST;
  typedef Tpetra::Map<LO, GO, NT> map_type;
  typedef Tpetra::Export<LO, GO, NT> export_type;
  typedef Tpetra::CrsMatrix<ST, LO, GO, NT> matrix_type;

  RCP<const Teuchos::Comm<int> > comm = rowMap->getComm ();
  const int myRank = comm->getRank ();

  const GST globalNumElts = rowMap->getGlobalNumElements ();
  const size_t myNumElts = (myRank == 0) ? as<size_t> (globalNumElts) : 0;
  RCP<const map_type> gatherRowMap = computeGatherMap (rowMap, rcpFromRef (out), dbg);
  matrix_type A_gather (gatherRowMap, 3);

  if (myRank == 0) {
    Array<GO> ind (3);
    Array<ST> val (3);
    for (size_t myRow = 0; myRow < myNumElts; ++myRow) {
      const GO globalRow = gatherRowMap->getGlobalElement (myRow);
      if (globalRow == gatherRowMap->getMinAllGlobalIndex ()) {
        val[0] = as<ST> (2);
        val[1] = as<ST> (-1);
        ind[0] = globalRow;
        ind[1] = globalRow + 1;
        A_gather.insertGlobalValues (globalRow, ind.view (0, 2), val.view (0, 2));
      }
      else if (globalRow == gatherRowMap->getMaxAllGlobalIndex ()) {
        val[0] = as<ST> (-1);
        val[1] = as<ST> (2);
        ind[0] = globalRow - 1;
        ind[1] = globalRow;
        A_gather.insertGlobalValues (globalRow, ind.view (0, 2), val.view (0, 2));
      }
      else {
        val[0] = as<ST> (-1);
        val[1] = as<ST> (2);
        val[2] = as<ST> (-1);
        ind[0] = globalRow - 1;
        ind[1] = globalRow;
        ind[2] = globalRow + 1;
        A_gather.insertGlobalValues (globalRow, ind.view (0, 3), val.view (0, 3));
      }
    }
  }
  A_gather.fillComplete (rowMap, rowMap);
  RCP<matrix_type> A = rcp (new matrix_type (rowMap, as<size_t> (0)));
  export_type exp (gatherRowMap, rowMap);
  A->doExport (A_gather, exp, Tpetra::INSERT);
  A->fillComplete (rowMap, rowMap);
  return A;
}

template<class ScalarType, class LocalOrdinalType, class GlobalOrdinalType, class NodeType>
bool
testCrsMatrix (Teuchos::FancyOStream& out, const GlobalOrdinalType indexBase)
{
  typedef ScalarType ST;
  typedef LocalOrdinalType LO;
  typedef GlobalOrdinalType GO;
  typedef NodeType NT;
  typedef Tpetra::Map<LO, GO, NT> map_type;
  typedef Tpetra::CrsMatrix<ST, LO, GO, NT> crs_matrix_type;
  bool result = true; // current Boolean result; reused below
  bool success = true;

  out << "Test: CrsMatrix Matrix Market I/O, w/ Map with index base "
      << indexBase << endl;
  OSTab tab1 (out);

  RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();

  out << "Original sparse matrix:" << endl;
  out << matrix_symRealSmall << endl;

  out << "Creating the row Map" << endl;
  const global_size_t globalNumElts = 5;
  RCP<const map_type> rowMap =
    rcp (new map_type (globalNumElts, indexBase, comm,
                       Tpetra::GloballyDistributed));

  out << "Reading in the matrix" << endl;
  std::istringstream inStr (matrix_symRealSmall);
  RCP<const map_type> colMap;
  RCP<const map_type> domainMap = rowMap;
  RCP<const map_type> rangeMap = rowMap;
  typedef Tpetra::MatrixMarket::Reader<crs_matrix_type> reader_type;
  RCP<crs_matrix_type> A =
    reader_type::readSparse (inStr, rowMap, colMap, domainMap, rangeMap,
                             callFillComplete, tolerant, debug);

  out << "Creating original matrix" << endl;
  RCP<crs_matrix_type> A_orig =
    createSymRealSmall<ST, LO, GO, NT> (rowMap, out, debug);

  out << "Comparing read-in matrix to original matrix" << endl;
  result = compareCrsMatrix<crs_matrix_type> (*A_orig, *A, out);
  bool local_success = true;
  TEUCHOS_TEST_EQUALITY( result, true, out, local_success );
  if (! result) { // see if ignoring zero values helps
    result = compareCrsMatrixValues<crs_matrix_type> (*A_orig, *A, out);
    local_success = true;
    TEUCHOS_TEST_EQUALITY( result, true, out, local_success );
  }
  success = success && local_success;

  out << "Writing out the original matrix" << endl;
  std::ostringstream outStr;
  typedef Tpetra::MatrixMarket::Writer<crs_matrix_type> writer_type;
  writer_type::writeSparse (outStr, A_orig, debug);

  out << "Reading it in again and comparing with original matrix" << endl;
  std::istringstream inStr2 (outStr.str ());
  RCP<crs_matrix_type> A_orig2 =
    reader_type::readSparse (inStr2, rowMap, colMap, domainMap, rangeMap,
                             callFillComplete, tolerant, debug);
  result = compareCrsMatrix<crs_matrix_type> (*A_orig, *A_orig2, out);
  local_success = true;
  TEUCHOS_TEST_EQUALITY( result, true, out, local_success );
  if (! result) { // see if ignoring zero values helps
    result = compareCrsMatrixValues<crs_matrix_type> (*A_orig, *A_orig2, out);
    local_success = true;
    TEUCHOS_TEST_EQUALITY( result, true, out, local_success );
  }
  success = success && local_success;

  return success;
}


template<class ScalarType, class LocalOrdinalType, class GlobalOrdinalType, class NodeType>
bool
testCrsMatrixPerFile (Teuchos::FancyOStream& out, const GlobalOrdinalType indexBase)
{
  typedef ScalarType ST;
  typedef LocalOrdinalType LO;
  typedef GlobalOrdinalType GO;
  typedef NodeType NT;
  typedef Tpetra::Map<LO, GO, NT> map_type;
  typedef Tpetra::CrsMatrix<ST, LO, GO, NT> crs_matrix_type;
  bool result = true; // current Boolean result; reused below
  bool success = true;

  out << "Test: CrsMatrix Per File Matrix Market I/O, w/ Map with index base "
      << indexBase << endl;
  OSTab tab1 (out);

  RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();

  out << "Creating the row Map" << endl;
  const global_size_t globalNumElts = 5;
  RCP<const map_type> rowMap =
    rcp (new map_type (globalNumElts, indexBase, comm,
                       Tpetra::GloballyDistributed));

  out << "Creating original matrix" << endl;
  RCP<crs_matrix_type> A_orig =
    createSymRealSmall<ST, LO, GO, NT> (rowMap, out, debug);

  out << "Original sparse matrix:" << endl;
  A_orig->describe(out, Teuchos::VERB_EXTREME);



  // We'll add the number of total MPI ranks to the suffix
  // so if ctest runs multiple jobs w/ different num procs at 
  // the same time, they won't clash   
  std::string prefix = "outfile_";
  std::string suffix = std::string("_") + std::to_string(comm->getSize());

  // Write out the matrix, rank by rank
  typedef Tpetra::MatrixMarket::Writer<crs_matrix_type> writer_type;
  writer_type::writeSparsePerRank(prefix,suffix,*A_orig,"Original","");

  // Read matrix back in
  typedef Tpetra::MatrixMarket::Reader<crs_matrix_type> reader_type;
  RCP<const map_type> row_map = rowMap;
  RCP<const map_type> col_map, domain_map = row_map, range_map = row_map;
  RCP<crs_matrix_type> A_new = reader_type::readSparsePerRank(prefix,suffix,row_map,col_map,domain_map,range_map);

  out << "New sparse matrix:" << endl;
  A_new->describe(out, Teuchos::VERB_EXTREME);

  out << "Comparing read-in matrix to original matrix" << endl;
  result = compareCrsMatrix<crs_matrix_type> (*A_orig, *A_new, out);
  bool local_success = true;
  TEUCHOS_TEST_EQUALITY( result, true, out, local_success );
  if (! result) { // see if ignoring zero values helps
    result = compareCrsMatrixValues<crs_matrix_type> (*A_orig, *A_new, out);
    local_success = true;
    TEUCHOS_TEST_EQUALITY( result, true, out, local_success );
  }
  success = success && local_success;

  return success;
}

template<class ScalarType, class LocalOrdinalType, class GlobalOrdinalType, class NodeType>
bool
testCrsMatrixPerFileNonContiguous (Teuchos::FancyOStream& out, const GlobalOrdinalType indexBase)
{
  typedef ScalarType ST;
  typedef LocalOrdinalType LO;
  typedef GlobalOrdinalType GO;
  typedef NodeType NT;
  typedef Tpetra::Map<LO, GO, NT> map_type;
  typedef Tpetra::CrsMatrix<ST, LO, GO, NT> crs_matrix_type;
  bool result = true; // current Boolean result; reused below
  bool success = true;

  out << "Test: Noncontiguous row map CrsMatrix Per File Matrix Market I/O, w/ Map with index base "
      << indexBase << endl;
  OSTab tab1 (out);

  RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();

  // Two rows per rank, reverse numbered
  out << "Creating the non-contiguous row Map" << endl;
  GO INVALID = Teuchos::OrdinalTraits<GO>::invalid();
  Teuchos::Array<GO> my_gids(2);
  my_gids[0] =  2*comm->getSize() - 2*comm->getRank() - 2 + indexBase;
  my_gids[1] = my_gids[0] + 1;
  RCP<const map_type> rowMap =
    rcp (new map_type (INVALID,my_gids(),indexBase,comm));

  out << "Creating original matrix" << endl;
  RCP<crs_matrix_type> A_orig =
    createSymRealSmall<ST, LO, GO, NT> (rowMap, out, debug);

  out << "Original sparse matrix:" << endl;
  A_orig->describe(out, Teuchos::VERB_EXTREME);

  // We'll add the number of total MPI ranks to the suffix
  // so if ctest runs multiple jobs w/ different num procs at 
  // the same time, they won't clash   
  std::string prefix = "outfile_";
  std::string suffix = std::string("_") + std::to_string(comm->getSize());

  // Write out the matrix, rank by rank
  typedef Tpetra::MatrixMarket::Writer<crs_matrix_type> writer_type;
  writer_type::writeSparsePerRank(prefix,suffix,*A_orig,"Original","");

  // Read matrix back in
  typedef Tpetra::MatrixMarket::Reader<crs_matrix_type> reader_type;
  RCP<const map_type> row_map = rowMap;
  RCP<const map_type> col_map, domain_map = row_map, range_map = row_map;
  RCP<crs_matrix_type> A_new = reader_type::readSparsePerRank(prefix,suffix,row_map,col_map,domain_map,range_map);

  out << "New sparse matrix:" << endl;
  A_new->describe(out, Teuchos::VERB_EXTREME);

  out << "Comparing read-in matrix to original matrix" << endl;
  result = compareCrsMatrix<crs_matrix_type> (*A_orig, *A_new, out);
  bool local_success = true;
  TEUCHOS_TEST_EQUALITY( result, true, out, local_success );
  if (! result) { // see if ignoring zero values helps
    result = compareCrsMatrixValues<crs_matrix_type> (*A_orig, *A_new, out);
    local_success = true;
    TEUCHOS_TEST_EQUALITY( result, true, out, local_success );
  }
  success = success && local_success;

  return success;
}



} // namespace (anonymous)

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrixOutputInput, IndexBase0, ST, LO, GO, NT )
{
  const GO indexBase = 0;
  success = testCrsMatrix<ST, LO, GO, NT> (out, indexBase);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrixOutputInput, IndexBase1, ST, LO, GO, NT )
{
  const GO indexBase = 1;
  success = testCrsMatrix<ST, LO, GO, NT> (out, indexBase);
}


TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrixPerFile, IndexBase0, ST, LO, GO, NT )
{
  const GO indexBase = 0;
  success = testCrsMatrixPerFile<ST, LO, GO, NT> (out, indexBase);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrixPerFile, IndexBase1, ST, LO, GO, NT )
{
  const GO indexBase = 1;
  success = testCrsMatrixPerFile<ST, LO, GO, NT> (out, indexBase);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrixPerFileNonContiguous, IndexBase0, ST, LO, GO, NT )
{
  const GO indexBase = 0;
  success = testCrsMatrixPerFileNonContiguous<ST, LO, GO, NT> (out, indexBase);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrixPerFileNonContiguous, IndexBase1, ST, LO, GO, NT )
{
  const GO indexBase = 1;
  success = testCrsMatrixPerFileNonContiguous<ST, LO, GO, NT> (out, indexBase);
}


// We instantiate tests for all combinations of the following parameters:
//   - indexBase = {0, 1}
//   - Scalar = {double, float}

#if defined(HAVE_TPETRA_INST_DOUBLE)
#  define UNIT_TEST_GROUP( LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrixOutputInput, IndexBase0, double, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrixOutputInput, IndexBase1, double, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrixPerFile, IndexBase0, double, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrixPerFile, IndexBase1, double, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrixPerFileNonContiguous, IndexBase0, double, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrixPerFileNonContiguous, IndexBase1, double, LO, GO, NODE ) 

#elif defined(HAVE_TPETRA_INST_FLOAT)
#  define UNIT_TEST_GROUP( LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrixOutputInput, IndexBase0, float, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrixOutputInput, IndexBase1, float, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrixPerFile, IndexBase0, float, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrixPerFile, IndexBase1, float, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrixPerFileNonContiguous, IndexBase0, float, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrixPerFileNonContiguous, IndexBase1, float, LO, GO, NODE ) 
#else
#  define UNIT_TEST_GROUP( LO, GO, NODE )
#endif


  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_LGN( UNIT_TEST_GROUP )

