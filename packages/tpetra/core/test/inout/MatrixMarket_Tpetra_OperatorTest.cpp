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

#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_Util.hpp> // sort2, merge2
#include <Teuchos_UnitTestHarness.hpp>
#include "TpetraCore_ETIHelperMacros.h"

namespace { // anonymous

using Teuchos::Array;
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
  typedef Tpetra::global_size_t GST;
  using Teuchos::arcp;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
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
    // that's what MPI uses.
    {
      // Make sure that the conversion from size_t to int won't
      // overflow on any process.  It's OK to do the all-reduce,
      // because this is a test, not performance-oriented code.
      const int lclSuccess =
        (oneToOneMap->getNodeNumElements () <= static_cast<size_t> (INT_MAX)) ?
        1 : 0;
      int gblSuccess = 1;
      using Teuchos::outArg;
      reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
      TEUCHOS_TEST_FOR_EXCEPTION
        (gblSuccess != 1, std::runtime_error, "On some process, "
         "oneToOneMap->getNodeNumElements() = "
         << oneToOneMap->getNodeNumElements () << " > INT_MAX = " << INT_MAX
         << ".  This means that the code below that uses MPI_Gather to collect "
         "indices onto Process 0 will be incorrect.");
    }

    const int myEltCount = static_cast<int> (oneToOneMap->getNodeNumElements ());
    Array<int> recvCounts (numProcs);
    const int rootProc = 0;
    gather<int, int> (&myEltCount, 1, recvCounts.getRawPtr (), 1, rootProc, *comm);

    ArrayView<const GO> myGlobalElts = oneToOneMap->getNodeElementList ();
    const int numMyGlobalElts = static_cast<int> (myGlobalElts.size ());
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
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    gatherMap =
      rcp (new MapType (INVALID, allElts, oneToOneMap->getIndexBase (), comm));
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
  typedef typename CrsMatrixType::scalar_type ST;
  typedef typename CrsMatrixType::global_ordinal_type GO;
  typedef typename ArrayView<const GO>::size_type size_type;

  Teuchos::OSTab tab (Teuchos::rcpFromRef (out));
  int localEqual = 1;

  //
  // Are my local matrices equal?
  //
  Array<GO> indOrig, ind;
  Array<ST> valOrig, val;
  size_t numEntriesOrig = 0;
  size_t numEntries = 0;

  ArrayView<const GO> localElts = A.getRowMap ()->getNodeElementList ();
  const size_type numLocalElts = localElts.size ();
  for (size_type i = 0; i < numLocalElts; ++i) {
    const GO globalRow = localElts[i];
    numEntriesOrig = A_orig.getNumEntriesInGlobalRow (globalRow);
    numEntries = A.getNumEntriesInGlobalRow (globalRow);

    if (numEntriesOrig != numEntries) {
      localEqual = 0;
      break;
    }
    indOrig.resize (numEntriesOrig);
    valOrig.resize (numEntriesOrig);
    A_orig.getGlobalRowCopy (globalRow, indOrig (), valOrig (), numEntriesOrig);
    ind.resize (numEntries);
    val.resize (numEntries);
    A.getGlobalRowCopy (globalRow, ind (), val (), numEntries);

    // Global row entries are not necessarily sorted.  Sort them so
    // we can compare them.
    Tpetra::sort2 (indOrig.begin (), indOrig.end (), valOrig.begin ());
    Tpetra::sort2 (ind.begin (), ind.end (), val.begin ());

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

  Teuchos::OSTab tab (Teuchos::rcpFromRef (out));

  //
  // Are my local matrices equal?
  //
  Array<GO> indOrig, ind;
  Array<ST> valOrig, val;
  size_t numEntriesOrig = 0;
  size_t numEntries = 0;

  ArrayView<const GO> localElts = A.getRowMap ()->getNodeElementList ();
  const size_type numLocalElts = localElts.size ();
  MT localDiff = STM::zero (); // \sum_{i,j} |A_orig(i,j) - A(i,j)| locally
  for (size_type i = 0; i < numLocalElts; ++i) {
    const GO globalRow = localElts[i];
    numEntriesOrig = A_orig.getNumEntriesInGlobalRow (globalRow);
    numEntries = A.getNumEntriesInGlobalRow (globalRow);

    indOrig.resize (numEntriesOrig);
    valOrig.resize (numEntriesOrig);
    A_orig.getGlobalRowCopy (globalRow, indOrig (), valOrig (), numEntriesOrig);
    ind.resize (numEntries);
    val.resize (numEntries);
    A.getGlobalRowCopy (globalRow, ind (), val (), numEntries);

    // Global row entries are not necessarily sorted.  Sort them
    // (and their values with them) so we can merge their values.
    Tpetra::sort2 (indOrig.begin (), indOrig.end (), valOrig.begin ());
    Tpetra::sort2 (ind.begin (), ind.end (), val.begin ());

    //
    // Merge repeated values in each set of indices and values.
    //

    typename Array<GO>::iterator indOrigIter = indOrig.begin ();
    typename Array<ST>::iterator valOrigIter = valOrig.begin ();
    typename Array<GO>::iterator indOrigEnd = indOrig.end ();
    typename Array<ST>::iterator valOrigEnd = valOrig.end ();
    Tpetra::merge2 (indOrigEnd, valOrigEnd, indOrigIter, indOrigEnd, valOrigIter, valOrigEnd);

    typename Array<GO>::iterator indIter = ind.begin ();
    typename Array<ST>::iterator valIter = val.begin ();
    typename Array<GO>::iterator indEnd = ind.end ();
    typename Array<ST>::iterator valEnd = val.end ();
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
  const size_t myNumElts = (myRank == 0) ?
    static_cast<size_t> (globalNumElts) : static_cast<size_t> (0);
  RCP<const map_type> gatherRowMap = computeGatherMap (rowMap, rcpFromRef (out), dbg);
  matrix_type A_gather (gatherRowMap, static_cast<size_t> (0));

  if (myRank == 0) {
    Array<GO> ind (3);
    Array<ST> val (3);

    const ST ONE = Teuchos::ScalarTraits<ST>::one ();
    const ST TWO = ONE + ONE;

    for (size_t myRow = 0; myRow < myNumElts; ++myRow) {
      const GO globalRow = gatherRowMap->getGlobalElement (myRow);
      if (globalRow == gatherRowMap->getMinAllGlobalIndex ()) {
        val[0] = TWO;
        val[1] = -ONE;
        ind[0] = globalRow;
        ind[1] = globalRow + 1;
        A_gather.insertGlobalValues (globalRow, ind.view (0, 2), val.view (0, 2));
      }
      else if (globalRow == gatherRowMap->getMaxAllGlobalIndex ()) {
        val[0] = -ONE;
        val[1] = TWO;
        ind[0] = globalRow - 1;
        ind[1] = globalRow;
        A_gather.insertGlobalValues (globalRow, ind.view (0, 2), val.view (0, 2));
      }
      else {
        val[0] = -ONE;
        val[1] = TWO;
        val[2] = -ONE;
        ind[0] = globalRow - 1;
        ind[1] = globalRow;
        ind[2] = globalRow + 1;
        A_gather.insertGlobalValues (globalRow, ind.view (0, 3), val.view (0, 3));
      }
    }
  }
  A_gather.fillComplete (rowMap, rowMap);
  RCP<matrix_type> A = rcp (new matrix_type (rowMap, static_cast<size_t> (0)));
  export_type exp (gatherRowMap, rowMap);
  A->doExport (A_gather, exp, Tpetra::INSERT);
  A->fillComplete (rowMap, rowMap);
  return A;
}

template<class ScalarType, class LocalOrdinalType, class GlobalOrdinalType, class NodeType>
bool
testCrsMatrix (Teuchos::FancyOStream& out, const GlobalOrdinalType indexBase)
{
  typedef Tpetra::global_size_t GST;
  typedef ScalarType ST;
  typedef LocalOrdinalType LO;
  typedef GlobalOrdinalType GO;
  typedef NodeType NT;
  typedef Tpetra::Map<LO, GO, NT> map_type;
  typedef Tpetra::CrsMatrix<ST, LO, GO, NT> crs_matrix_type;
  bool result = true; // current Boolean result; reused below
  bool success = true;

  out << "Test: Operator Matrix Market I/O, w/ Map with index base "
      << indexBase << endl;
  OSTab tab1 (out);

  RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();

  out << "Creating the row Map" << endl;
  const GST globalNumElts = 34;
  RCP<const map_type> rowMap =
    rcp (new map_type (globalNumElts, indexBase, comm,
                       Tpetra::GloballyDistributed));

  RCP<const map_type> domainMap = rowMap;
  RCP<const map_type> rangeMap = rowMap;
  typedef Tpetra::MatrixMarket::Reader<crs_matrix_type> reader_type;

  out << "Creating original matrix" << endl;
  RCP<crs_matrix_type> A_orig =
    createSymRealSmall<ST, LO, GO, NT> (rowMap, out, debug);

  typedef Tpetra::MatrixMarket::Writer<crs_matrix_type> op_writer_type;
  Teuchos::ParameterList pl;
  pl.set("probing size",5);
  //pl.set("precision",12);
  //pl.set("zero-based indexing",true);
  //pl.set("print MatrixMarket header",false);

  // mfh 23 May 2015: It's unwise to write to files in tests.
  // Instead, we write to an output stream.

  out << "Writing out the original Operator to an output stream" << endl;
  std::ostringstream originalOperatorFile;
  op_writer_type::writeOperator (originalOperatorFile, *A_orig, pl);

  out << "Reading it in again and comparing with original matrix" << endl;
  std::istringstream readInMatrixFile (originalOperatorFile.str ());

  RCP<const map_type> colMap;
  RCP<crs_matrix_type> A_orig2 =
    reader_type::readSparse (readInMatrixFile, rowMap, colMap,
                             domainMap, rangeMap,
                             callFillComplete, tolerant, debug);
  result = compareCrsMatrix<crs_matrix_type> (*A_orig, *A_orig2, out);
  bool local_success = true;
  TEUCHOS_TEST_EQUALITY( result, true, out, local_success );
  if (! result) { // see if ignoring zero values helps
    result = compareCrsMatrixValues<crs_matrix_type> (*A_orig, *A_orig2, out);
    local_success = true;
    TEUCHOS_TEST_EQUALITY( result, true, out, local_success );
  }
  success = success && local_success;

  return success;
}

} // namespace (anonymous)

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( OperatorOutput, IndexBase0, ST, LO, GO, NT )
{
  const GO indexBase = 0;
  success = testCrsMatrix<ST, LO, GO, NT> (out, indexBase);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( OperatorOutput, IndexBase1, ST, LO, GO, NT )
{
  const GO indexBase = 1;
  success = testCrsMatrix<ST, LO, GO, NT> (out, indexBase);
}

// We instantiate tests for all combinations of the following parameters:
//   - indexBase = {0, 1}
//   - Scalar = {double, float}

#if defined(HAVE_TPETRA_INST_DOUBLE)
#  define UNIT_TEST_GROUP( LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( OperatorOutput, IndexBase0, double, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( OperatorOutput, IndexBase1, double, LO, GO, NODE )

#elif defined(HAVE_TPETRA_INST_FLOAT)
#  define UNIT_TEST_GROUP( LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( OperatorOutput, IndexBase0, float, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( OperatorOutput, IndexBase1, float, LO, GO, NODE )

#else
#  define UNIT_TEST_GROUP( LO, GO, NODE )
#endif


  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_LGN( UNIT_TEST_GROUP )

