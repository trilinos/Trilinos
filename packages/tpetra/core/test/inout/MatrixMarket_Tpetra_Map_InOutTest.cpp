// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <TpetraCore_ETIHelperMacros.h>

#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_Core.hpp>
#include <Teuchos_UnitTestHarness.hpp>

namespace { // anonymous

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
using std::cerr;
using std::endl;

typedef Tpetra::global_size_t GST;

// Whether to print copious debugging output to stderr when doing
// Matrix Market input and output.  This affects all the tests.
const bool debug = false;


TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MapOutputInput, ContigUniformIndexBase0, LO, GO )
{
  typedef Tpetra::Map<LO, GO> map_type;
  const bool tolerant = false;
  const bool globallyVerbose = true;

  RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();

  const size_t localNumElts = 10;
  const GST globalNumElts = localNumElts * static_cast<GST> (numProcs);
  const GO indexBase = 0;

  if (myRank == 0) {
    std::ostringstream os;
    os << "Test: Output a Tpetra::Map to a MatrixMarket file" << endl
       << "  - Contiguous, uniform, index base " << indexBase << endl
       << "  - Index base 0" << endl
       << "  - " << numProcs << " process(es)" << endl;
    cerr << os.str ();
  }

  if (globallyVerbose) {
    std::ostringstream os;
    os << "Proc " << myRank << ": Creating the Map" << endl;
    cerr << os.str ();
  }
  map_type map (globalNumElts, indexBase, comm, Tpetra::GloballyDistributed);

  // Write to a distinct output stream, so we can read it back in again.
  if (globallyVerbose) {
    std::ostringstream os;
    os << "Proc " << myRank << ": Writing Map to output stream" << endl;
    cerr << os.str ();
  }
  std::ostringstream mapOutStream;
  // The Scalar type doesn't matter, since we're just writing a Map.
  typedef Tpetra::CrsMatrix<double, LO, GO> crs_matrix_type;
  typedef Tpetra::MatrixMarket::Writer<crs_matrix_type> writer_type;
  writer_type::writeMap (mapOutStream, map, false && globallyVerbose);

  if (globallyVerbose) {
    std::ostringstream os;
    os << "Proc " << myRank << ": Done writing Map" << endl;
    cerr << os.str ();
  }

  // Make sure that all processes finished the above step before
  // continuing the test.  This ensures that if anything above is
  // broken and has "dangling messages," they won't get mixed up with
  // readMap() below.
  comm->barrier ();

  if (myRank == 0) {
    std::ostringstream os;
    os << "Result of writing the Map:" << endl
       << mapOutStream.str () << endl;
    cerr << os.str ();
  }

  comm->barrier ();

  //
  // Read the Map back in from the output stream to which it was written.
  //

  if (globallyVerbose) {
    std::ostringstream os;
    os << "Proc " << myRank << ": Reading Map back in" << endl;
    cerr << os.str ();
  }
  // This input stream can exist on all processes,
  // but only Process 0 will read from it.
  std::istringstream mapInStream (mapOutStream.str ());
  typedef Tpetra::MatrixMarket::Reader<crs_matrix_type> reader_type;
  RCP<const map_type> inMap =
    reader_type::readMap (mapInStream, map.getComm (),
                          tolerant, globallyVerbose);
  if (globallyVerbose) {
    std::ostringstream os;
    os << "Proc " << myRank << ": Done reading Map" << endl;
    cerr << os.str ();
  }

  // Make sure that all processes finished the above step before
  // continuing the test.  This ensures that if anything above is
  // broken and has "dangling messages," they won't get mixed up with
  // isSameAs() below.
  comm->barrier ();

  if (globallyVerbose) {
    std::ostringstream os;
    os << "Proc " << myRank << ": Testing sameness" << endl;
    cerr << os.str ();
  }
  TEST_ASSERT( map.isSameAs (*inMap) );

  // Now try readMap again.

  // This input stream can exist on all processes,
  // but only Process 0 will read from it.
  std::istringstream mapInStream2 (mapOutStream.str ());
  RCP<const map_type> inMap2 =
    reader_type::readMap (mapInStream2, map.getComm (),
                          tolerant, globallyVerbose);
  // Make sure that all processes finished the above step before
  // continuing the test.
  comm->barrier ();

  if (globallyVerbose) {
    std::ostringstream os;
    os << "Proc " << myRank << ": Testing sameness" << endl;
    cerr << os.str ();
  }
  TEST_ASSERT( map.isSameAs (*inMap2) );

  if (globallyVerbose) {
    std::ostringstream os;
    os << "Proc " << myRank << ": Finished with test" << endl;
    cerr << os.str ();
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MapOutputInput, ContigUniformIndexBase1, LO, GO )
{
  typedef Tpetra::Map<LO, GO> map_type;

  const bool tolerant = false;
  out << "Test: output a contiguous uniform Tpetra::Map (index base 1) "
    "to a Matrix Market file" << endl;
  OSTab tab1 (out);

  RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();
  // Forestall compiler warnings about unused variables.
  (void) myRank;
  (void) numProcs;

  const size_t localNumElts = 10;
  const GST globalNumElts = localNumElts * static_cast<GST> (numProcs);
  const GO indexBase = 1;

  out << "Creating the Map" << endl;
  map_type map (globalNumElts, indexBase, comm, Tpetra::GloballyDistributed);

  // Write to a distinct output stream, so we can read it back in again.
  out << "Writing Map to output stream" << endl;
  std::ostringstream mapOutStream;
  // The Scalar type doesn't matter, since we're just writing a Map.
  typedef Tpetra::CrsMatrix<double, LO, GO> crs_matrix_type;
  typedef Tpetra::MatrixMarket::Writer<crs_matrix_type> writer_type;
  writer_type::writeMap (mapOutStream, map, debug);
  out << "Result of writing the Map:" << endl;
  if (myRank == 0) {
    OSTab tab2 (out);
    out << mapOutStream.str () << endl;
  }

  out << "Reading Map back in from the output stream" << endl;
  std::istringstream mapInStream (mapOutStream.str ());
  typedef Tpetra::MatrixMarket::Reader<crs_matrix_type> reader_type;
  RCP<const map_type> inMap =
    reader_type::readMap (mapInStream, map.getComm (),
                          tolerant, debug);

  TEST_EQUALITY(map.isSameAs (*inMap), true);
}


TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MapOutputInput, ContigNonuniformIndexBase0, LO, GO )
{
  typedef Tpetra::Map<LO, GO> map_type;

  const bool tolerant = false;
  out << "Test: output a contiguous nonuniform Tpetra::Map (index base 0) "
    "to a Matrix Market file" << endl;
  OSTab tab1 (out);

  RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
  const int myRank = comm->getRank();
  const int numProcs = comm->getSize();

  // Proc 0 gets a different number of local elements, so that the Map is nonuniform.
  const size_t localNumElts = (myRank == 0) ? 11 : 10;
  const GST globalNumElts = 10 * static_cast<GST> (numProcs) + 1;
  const GO indexBase = 0;

  out << "Creating the Map" << endl;
  map_type map (globalNumElts, localNumElts, indexBase, comm);

  // Write to a distinct output stream, so we can read it back in again.
  out << "Writing Map to output stream" << endl;
  std::ostringstream mapOutStream;
  // The Scalar type doesn't matter, since we're just writing a Map.
  typedef Tpetra::CrsMatrix<double, LO, GO> crs_matrix_type;
  typedef Tpetra::MatrixMarket::Writer<crs_matrix_type> writer_type;
  writer_type::writeMap (mapOutStream, map, debug);
  out << "Result of writing the Map:" << endl;
  if (myRank == 0) {
    OSTab tab2 (out);
    out << mapOutStream.str () << endl;
  }

  out << "Reading Map back in from the output stream" << endl;
  std::istringstream mapInStream (mapOutStream.str ());
  typedef Tpetra::MatrixMarket::Reader<crs_matrix_type> reader_type;
  RCP<const map_type> inMap =
    reader_type::readMap (mapInStream, map.getComm (), 
                          tolerant, debug);

  TEST_EQUALITY(map.isSameAs (*inMap), true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MapOutputInput, ContigNonuniformIndexBase1, LO, GO )
{
  typedef Tpetra::Map<LO, GO> map_type;

  const bool tolerant = false;
  out << "Test: output a contiguous nonuniform Tpetra::Map (index base 1) "
    "to a Matrix Market file" << endl;
  OSTab tab1 (out);

  RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();

  // Proc 0 gets a different number of local elements, so that the Map is nonuniform.
  const size_t localNumElts = (myRank == 0) ? 11 : 10;
  const GST globalNumElts = 10 * static_cast<GST> (numProcs) + 1;
  const GO indexBase = 1;

  out << "Creating the Map" << endl;
  map_type map (globalNumElts, localNumElts, indexBase, comm);

  // Write to a distinct output stream, so we can read it back in again.
  out << "Writing Map to output stream" << endl;
  std::ostringstream mapOutStream;
  // The Scalar type doesn't matter, since we're just writing a Map.
  typedef Tpetra::CrsMatrix<double, LO, GO> crs_matrix_type;
  typedef Tpetra::MatrixMarket::Writer<crs_matrix_type> writer_type;
  writer_type::writeMap (mapOutStream, map, debug);
  out << "Result of writing the Map:" << endl;
  if (myRank == 0) {
    OSTab tab2 (out);
    out << mapOutStream.str () << endl;
  }

  out << "Reading Map back in from the output stream" << endl;
  std::istringstream mapInStream (mapOutStream.str ());
  typedef Tpetra::MatrixMarket::Reader<crs_matrix_type> reader_type;
  RCP<const map_type> inMap =
    reader_type::readMap (mapInStream, map.getComm (), 
                          tolerant, debug);

  TEST_EQUALITY(map.isSameAs (*inMap), true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MapOutputInput, NoncontigIndexBase0, LO, GO )
{
  typedef Tpetra::Map<LO, GO> map_type;

  RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();

  const bool tolerant = false;
  out << "Test: output a noncontiguous Tpetra::Map (index base 0) "
    "to a Matrix Market file" << endl;
  OSTab tab1 (out);

  // Just to make sure that the Map is really stored noncontiguously,
  // we only have it contain even-numbered GIDs.
  const size_t localNumElts = 10;
  const GST globalNumElts = 10 * static_cast<GST> (numProcs);
  Array<GO> localGids (localNumElts);
  const GO indexBase = 0;

  const GO localStartGid = indexBase + 2 * as<GO> (myRank) * as<GO> (localNumElts);
  for (size_t k = 0; k < localNumElts; ++k) {
    localGids[k] = localStartGid + as<GO> (2*k);
  }
  out << "Creating the Map" << endl;
  map_type map (globalNumElts, localGids (), indexBase, comm);

  // Write to a distinct output stream, so we can read it back in again.
  out << "Writing Map to output stream" << endl;
  std::ostringstream mapOutStream;
  // The Scalar type doesn't matter, since we're just writing a Map.
  typedef Tpetra::CrsMatrix<double, LO, GO> crs_matrix_type;
  typedef Tpetra::MatrixMarket::Writer<crs_matrix_type> writer_type;
  writer_type::writeMap (mapOutStream, map, debug);
  out << "Result of writing the Map:" << endl;
  if (myRank == 0) {
    OSTab tab2 (out);
    out << mapOutStream.str () << endl;
  }

  out << "Reading Map back in from the output stream" << endl;
  std::istringstream mapInStream (mapOutStream.str ());
  typedef Tpetra::MatrixMarket::Reader<crs_matrix_type> reader_type;
  RCP<const map_type> inMap =
    reader_type::readMap (mapInStream, map.getComm (), 
                          tolerant, debug);

  TEST_EQUALITY(map.isSameAs (*inMap), true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MapOutputInput, NoncontigIndexBase1, LO, GO )
{
  typedef Tpetra::Map<LO, GO> map_type;

  RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();

  const bool tolerant = false;
  out << "Test: output a noncontiguous Tpetra::Map (index base 1) "
    "to a Matrix Market file" << endl;
  OSTab tab1 (out);

  // Just to make sure that the Map is really stored noncontiguously,
  // we only have it contain even-numbered GIDs.
  const size_t localNumElts = 10;
  const GST globalNumElts = 10 * static_cast<GST> (numProcs);
  Array<GO> localGids (localNumElts);
  const GO indexBase = 1;

  const GO localStartGid = indexBase + 2 * as<GO> (myRank) * as<GO> (localNumElts);
  for (size_t k = 0; k < localNumElts; ++k) {
    localGids[k] = localStartGid + as<GO> (2*k);
  }
  out << "Creating the Map" << endl;
  map_type map (globalNumElts, localGids (), indexBase, comm);

  // Write to a distinct output stream, so we can read it back in again.
  out << "Writing Map to output stream" << endl;
  std::ostringstream mapOutStream;
  // The Scalar type doesn't matter, since we're just writing a Map.
  typedef Tpetra::CrsMatrix<double, LO, GO> crs_matrix_type;
  typedef Tpetra::MatrixMarket::Writer<crs_matrix_type> writer_type;
  writer_type::writeMap (mapOutStream, map, debug);
  out << "Result of writing the Map:" << endl;
  if (myRank == 0) {
    OSTab tab2 (out);
    out << mapOutStream.str () << endl;
  }

  out << "Reading Map back in from the output stream" << endl;
  std::istringstream mapInStream (mapOutStream.str ());
  typedef Tpetra::MatrixMarket::Reader<crs_matrix_type> reader_type;
  RCP<const map_type> inMap =
    reader_type::readMap (mapInStream, map.getComm (), 
                          tolerant, debug);

  TEST_EQUALITY(map.isSameAs (*inMap), true);
}

// Noncontiguous, overlapping Map with index base 0.
TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( MapOutputInput, NoncontigOvrlpngIndBase0, LO, GO )
{
  typedef Tpetra::Map<LO, GO> map_type;
  RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();

  const bool tolerant = false;
  out << "Test: output a noncontiguous, overlapping Tpetra::Map "
    "(index base 0) to a Matrix Market file" << endl;
  OSTab tab1 (out);

  const size_t localNumElts = (myRank == 0 || myRank == numProcs - 1) ? 2 : 3;
  const GST globalNumElts = (numProcs <= 2) ?
    static_cast<GST> (numProcs * 2) :
    static_cast<GST> (4 + (numProcs - 2) * 3);
  Array<GO> localGids (localNumElts);
  const GO indexBase = 0;

  // Proc 0 has [0 1].
  // Proc 1 has [1 2 3].
  // Proc 2 has [2 3 4].
  // ...
  // Proc P-1 has [P-1 P].
  for (size_t k = 0; k < localNumElts; ++k) {
    localGids[k] = static_cast<GO> (myRank) + static_cast<GO> (k);
  }

  out << "Creating the Map" << endl;
  map_type map (globalNumElts, localGids (), indexBase, comm);

  // Write to a distinct output stream, so we can read it back in again.
  out << "Writing Map to output stream" << endl;
  std::ostringstream mapOutStream;
  // The Scalar type doesn't matter, since we're just writing a Map.
  typedef Tpetra::CrsMatrix<double, LO, GO> crs_matrix_type;
  typedef Tpetra::MatrixMarket::Writer<crs_matrix_type> writer_type;
  writer_type::writeMap (mapOutStream, map, debug);
  out << "Result of writing the Map:" << endl;
  if (myRank == 0) {
    OSTab tab2 (out);
    out << mapOutStream.str () << endl;
  }

  out << "Reading Map back in from the output stream" << endl;
  std::istringstream mapInStream (mapOutStream.str ());
  typedef Tpetra::MatrixMarket::Reader<crs_matrix_type> reader_type;
  RCP<const map_type> inMap =
    reader_type::readMap (mapInStream, map.getComm (), 
                          tolerant, debug);

  TEST_EQUALITY(map.isSameAs (*inMap), true);
}


// We instantiate tests for all combinations of the following parameters:
//   - All enabled (LO, GO) type combinations
//   - indexBase = {0, 1}
//   - {contiguous uniform, contiguous nonuniform, noncontiguous} Map

#define UNIT_TEST_GROUP( LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MapOutputInput, ContigUniformIndexBase0, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MapOutputInput, ContigUniformIndexBase1, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MapOutputInput, ContigNonuniformIndexBase0, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MapOutputInput, ContigNonuniformIndexBase1, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MapOutputInput, NoncontigIndexBase0, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MapOutputInput, NoncontigIndexBase1, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( MapOutputInput, NoncontigOvrlpngIndBase0, LO, GO )


TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_LG( UNIT_TEST_GROUP )

// mfh 05 Sep 2014: Must close namespace here for some reason, not
// above the instantiations.
} // namespace (anonymous)



