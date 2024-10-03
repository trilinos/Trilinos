// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_Util.hpp> // sort2
#include <Teuchos_UnitTestHarness.hpp>

using Tpetra::global_size_t;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::as;
using Teuchos::Comm;
using Teuchos::null;
using Teuchos::OSTab;
using Teuchos::ParameterList;
using Teuchos::ptr;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::rcpFromRef;
using Teuchos::REDUCE_MAX;
using Teuchos::REDUCE_MIN;
using Teuchos::reduceAll;
using std::endl;


namespace {

// Command-line arguments
std::string mapOutputFilename;
std::string graphOutputFilename;
std::string mapInputFilename;
std::string graphInputFilename;
bool debug = false;

TEUCHOS_STATIC_SETUP()
{
  Teuchos::CommandLineProcessor& clp = Teuchos::UnitTestRepository::getCLP ();
  clp.setOption ("mapOutputFilename", &mapOutputFilename,
                 "Name of file to which to write temporary Map output.  "
                 "If not provided, we will not test Map output.");
  clp.setOption ("graphOutputFilename", &graphOutputFilename,
                 "Name of file to which to write temporary CrsGraph output.  "
                 "If not provided, we write output to an std::ostringstream in "
                 "memory.  WARNING: If you don't provide an output file name, "
                 "the test may run out of memory if the input file is large, "
                 "since in that case we have to store two copies of the "
                 "graph's data in memory at once.");
  clp.setOption ("mapInputFilename", &mapInputFilename,
                 "Name of file from which to read the initial Map");
  clp.setOption ("graphInputFilename", &graphInputFilename,
                 "Name of file from which to read the initial CrsGraph "
                 "(using the input Map as its row Map)");
  clp.setOption ("debug", "release", &debug, "If true, print copious debugging "
                 "output to stderr on all processes.");
}


// Input graphs must be fill complete.
template<class CrsGraphType>
bool
compareCrsGraphMaps (const CrsGraphType& A_orig, const CrsGraphType& A, Teuchos::FancyOStream& out)
{
  Teuchos::OSTab tab (out);

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

// Input graphs must be fill complete, and all four of their Maps
// (row, column, domain, and range) must be the same.
template<class CrsGraphType>
bool
compareCrsGraph (const CrsGraphType& A_orig, const CrsGraphType& A, Teuchos::FancyOStream& out)
{
  typedef typename CrsGraphType::global_ordinal_type GO;
  typedef typename ArrayView<const GO>::size_type size_type;
  typedef typename CrsGraphType::nonconst_global_inds_host_view_type gids_type;

  OSTab tab (out);
  int localEqual = 1;

  //
  // Are my local graphs equal?
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
    Tpetra::sort(indOrig, indOrig.extent(0));
    Tpetra::sort(ind, ind.extent(0));

    for (size_t k = 0; k < numEntries; ++k) {
      // Values should be _exactly_ equal.
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



bool
testReadAndWriteFile (Teuchos::FancyOStream& out,
                      const std::string& mapOutFile,
                      const std::string& graphOutFile,
                      const std::string& mapInFile,
                      const std::string& graphInFile)
{
  using ST = double;
  using map_type = Tpetra::Map<>;
  using crs_matrix_type = Tpetra::CrsMatrix<ST>;
  using crs_graph_type = Tpetra::CrsGraph<>;
  using reader_type = Tpetra::MatrixMarket::Reader<crs_matrix_type>;
  using writer_type = Tpetra::MatrixMarket::Writer<crs_matrix_type>;
  const bool tolerant = false;

  bool result = true; // current Boolean result; reused below
  bool success = true; // used by TEST_EQUALITY

  out << "Test: Matrix Market, CrsGraph from row Map" << endl;
  OSTab tab1 (out);

  RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
  const int myRank = comm->getRank ();

  out << "Reading the (row) Map" << endl;
  RCP<const map_type> rowMap =
    reader_type::readMapFile (mapInFile, comm, tolerant, debug);

  if (mapOutFile != "") {
    out << "Writing the Map" << endl;
    writer_type::writeMapFile (mapOutFile, *rowMap);

    out << "Reading the written Map back in" << endl;
    RCP<const map_type> rowMap_out =
      reader_type::readMapFile (mapOutFile, comm, tolerant, debug);

    out << "Comparing the two Maps" << endl;
    result = rowMap->isSameAs (*rowMap_out);
    TEST_EQUALITY( result, true );
  }

  out << "Reading the CrsGraph" << endl;
  const bool callFillComplete = true;
  RCP<const map_type> colMap;
  RCP<const map_type> domainMap = rowMap;
  RCP<const map_type> rangeMap = rowMap;
  RCP<crs_graph_type> A_in =
    reader_type::readSparseGraphFile (graphInFile,
                                 rowMap, colMap, domainMap, rangeMap,
                                 callFillComplete, tolerant, debug);
  RCP<std::ostream> graphOutStream;
  if (graphOutFile == "") {
    if (myRank == 0) {
      out << "Process 0: No CrsGraph output file provided; storing output in "
        "memory." << endl << "WARNING: This may run out of memory if the input "
        "file is large, since we have to store two copies of the graph's data "
        "in memory at once." << endl;
      graphOutStream = rcp (new std::ostringstream);
    }
  } else {
    if (myRank == 0) {
      out << "Opening CrsGraph output file \"" << graphOutFile
          << "\" on Process 0 only" << endl;
      graphOutStream = rcp (new std::ofstream (graphOutFile.c_str ()));
    }
  }
  out << "Executing barrier to ensure no deadlock" << endl;
  comm->barrier ();
  out << "Writing the CrsGraph" << endl;
  writer_type::writeSparseGraph (*graphOutStream, *A_in, debug);

  out << "Reading the written CrsGraph back in" << endl;
  RCP<std::istream> graphInStream;
  if (graphOutFile == "") {
    if (myRank == 0) {
      // We have to store two copies of the output in memory, so that
      // we can construct the input stream from the output stream.
      RCP<std::ostringstream> matOss =
        rcp_dynamic_cast<std::ostringstream> (graphOutStream, true);
      graphInStream = rcp (new std::istringstream (matOss->str ()));
      graphOutStream = null;
    }
  } else {
    if (myRank == 0) {
      out << "Process 0: Closing CrsGraph output file" << endl;
      graphOutStream = null;
      out << "Process 0: Reopening CrsGraph output file \""
          << graphOutFile << "\"" << endl;
      graphInStream = rcp (new std::ifstream (graphOutFile.c_str ()));
    }
  }

  RCP<crs_graph_type> A_out =
    reader_type::readSparseGraph (*graphInStream,
                             rowMap, colMap, domainMap, rangeMap,
                             callFillComplete, tolerant, debug);

  out << "Test the two graphs for exact equality of structure" << endl;
  result = compareCrsGraph<crs_graph_type> (*A_in, *A_out, out);
  {
    OSTab tab (out);
    out << "- Graphs are " << (result ? "" : "NOT ") << "equal" << endl;
  }
  TEST_EQUALITY( result, true );

  // mfh 24 Feb 2014: Apparently, "set but not used" errors may still
  // show up on some compilers, even with use of the "(void) success"
  // idiom.  We can fix this by noting that success has to be true in
  // order for this function to return true, so we can just Boolean
  // AND the return value with success.
  //
  //(void) success; // silence compile warning ("set but not used")
  //return result;
  return result && success;
}

} // namespace (anonymous)


TEUCHOS_UNIT_TEST( MatrixMarket, CrsGraphFileTest )
{
  success = testReadAndWriteFile (out, mapOutputFilename, graphOutputFilename,
                                  mapInputFilename, graphInputFilename);
}



