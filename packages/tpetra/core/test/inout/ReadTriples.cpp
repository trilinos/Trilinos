// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_Details_ReadTriples.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include "Tpetra_Map.hpp"
#include "Teuchos_CommHelpers.hpp"
#ifdef HAVE_TPETRACORE_MPI
#  include "Teuchos_DefaultMpiComm.hpp"
#endif // HAVE_TPETRACORE_MPI
#include <complex>
#include <map>

namespace { // (anonymous)

using Teuchos::Comm;
using Teuchos::RCP;
using Teuchos::rcp;
using std::endl;
typedef ::Tpetra::Map<>::global_ordinal_type GO;

// Type of each (row index, column index) pair in the std::map we use
//   in this test for representing a sparse matrix.
template<class IndexType>
struct CooGraphEntry {
  IndexType row;
  IndexType col;
};

// Function comparing two CooGraphEntry structs, lexicographically,
//   first by row index, then by column index.  We use this as the
//   comparator in the std::map we use in this test for representing a
//   sparse matrix.
template<class IndexType>
struct CompareCooGraphEntries {
  bool
  operator () (const CooGraphEntry<IndexType>& x,
               const CooGraphEntry<IndexType>& y) const
  {
    return (x.row < y.row) || (x.row == y.row && x.col < y.col);
  }
};

// "Fake" file data to read in, for real-valued matrix entries.  I
// chose the indices so that we can tell if "symmetrization"
// (inserting both (i,j) and (j,i) when encountering (i,j)) works.
// The rule for computing the value corresponding to entry (i,j) is
// i*10 + j.
const char inputFileReal[] =
  "1 2 12.0\n"
  "3 4 34.0\n"
  "5 6 56.0\n"
  "7 8 78.0\n"
  "# this is a comment\n"
  "9 10 100.0\n"
  "11 12 122.0\n"
  "13 14 144.0\n";

// More "fake" file data to read in, for real-valued matrix entries.
// I picked this one so that on the last round, each receiving process
// (only 1 in this case: numProcs=2) has exactly
// maxNumEntriesPerMessage=2 entries to receive.  In the previous test
// case above, the last receiving process' last message has fewer
// entries than that to receive.  It's important to test both cases,
// so we make sure that the sending Process 0 doesn't send an extra
// message in the numEnt=0 case.
//
// I also put the comment line in the middle of a chunk, instead of
// between chunks.  (Chunk size is still 2 for this test case.)
const char inputFileReal2[] =
  "1 2 12.0\n"
  "3 4 34.0\n"
  "5 6 56.0\n"
  "# this is a comment\n"
  "7 8 78.0\n"
  "9 10 100.0\n"
  "11 12 122.0\n";

// Third real test case exercises a single chunk, that Process 0 can
// handle itself.  This makes sure that there aren't any extra
// messages that could cause hangs.  We also include a comment at the
// start of the "file," just to make sure that this works.
const char inputFileReal3[] =
  "# this is a comment\n"
  "1 2 12.0\n"
  "3 4 34.0\n";

// "Fake" file data to read in, for complex-valued matrix entries.  I
// chose the indices so that we can tell if "symmetrization"
// (inserting both (i,j) and (j,i) when encountering (i,j)) works, and
// chose the values so that we can tell if "Hermitian-ization" works.
//
// The rule for computing the real value corresponding to entry (i,j)
// is i*10 + j, and the complex value is 100 times the negative of
// that.
const char inputFileComplex[] =
  "1 2 12.0 -1200.0\n"
  "3 4 34.0 -3400.0\n"
  "5 6 56.0 -5600.0\n"
  "7 8 78.0 -7800.0\n"
  "# this is a comment\n"
  "9 10 100.0 -10000.0\n"
  "11 12 122.0 -12200.0\n"
  "13 14 144.0 -14400.0\n";

void
testReal (bool& success,
          Teuchos::FancyOStream& out,
          const char testCase[],
          const GO expectedRowInds[], // on the calling process
          const GO expectedColInds[], // on the calling process
          const double expectedVals[], // on the calling process
          const std::size_t expectedNumEnt, // on the calling process
          const Comm<int>& comm,
          const bool debug = false,
          std::ostream* errStrm = NULL)
{
  using Tpetra::Details::readAndDealOutTriples;
  using Teuchos::outArg;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  typedef double SC;
  // Type representing the local entries of a sparse matrix.  The set
  // sorts its entries lexicographically, first by row index, then by
  // column index.
  typedef std::map<CooGraphEntry<GO>, SC,
                   CompareCooGraphEntries<GO> > matrix_type;
  const int myRank = comm.getRank ();
  int lclSuccess = 1;
  int gblSuccess = 0; // output argument

  out << "Test readAndDealOutTriples for real-valued matrix entries" << endl;
  Teuchos::OSTab tab1 (out);

  // The entries of the sparse matrix that live on the calling process.
  matrix_type myMatrixEntries;
  // Closure that inserts a sparse matrix entry into myMatrixEntries
  // if it does not yet exist, else adds the input value to that of
  // the existing entry.
  std::function<int (const GO, const GO, const SC&)> insertEntry =
    [&myMatrixEntries] (const GO rowInd, const GO colInd, const SC& val) {
    auto iter = myMatrixEntries.find ({rowInd, colInd});
    if (iter == myMatrixEntries.end ()) {
      myMatrixEntries.insert ({{rowInd, colInd}, val});
    }
    else {
      iter->second += val;
    }
    return 0;
  };

  out << "Open the \"file\" on Process 0" << endl;
  std::istringstream inputStream (testCase);

  std::size_t curLineNum = 0;
  std::size_t totalNumEntriesRead = 0;

  // Keep this small, so that we can observe complete function
  // behavior.  If not less than the number of input lines, then
  // Process 0 will simply read all the entries and handle them
  // itself.  The tests below are hard-coded to assume this is 2, so
  // please change the tests if you want to change this number.
  const std::size_t maxNumEntriesPerMessage = 2;

  // One line of the input "file" is a comment line.  Note that these
  // values depend on the input "file."
  //const std::size_t expectedNumLinesRead = 8;
  //const std::size_t expectedTotalNumEntriesRead = 7;

  const bool tolerant = false;
  const int errCode =
    readAndDealOutTriples (inputStream, curLineNum, totalNumEntriesRead,
                           insertEntry, maxNumEntriesPerMessage, comm,
                           tolerant, errStrm, debug);
  TEST_EQUALITY( errCode, 0 );

  lclSuccess = success ? 1 : 0;
  gblSuccess = 0; // output argument
  reduceAll<int, int> (comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY( gblSuccess, 1 );
  if (gblSuccess != 1) {
    out << "Some process reported a nonzero error code!" << endl;
    return;
  }

  {
    std::ostringstream os;
    os << "Proc " << myRank << ": Got matrix entries: [";
    std::size_t k = 0;
    for (auto iter = myMatrixEntries.begin (); iter != myMatrixEntries.end (); ++iter, ++k) {
      os << "{row=" << iter->first.row << ", col=" << iter->first.col << ", val=" << iter->second << "}";
      if (k + 1 < myMatrixEntries.size ()) {
        os << ", ";
      }
    }
    os << "]" << endl;
    Tpetra::Details::gathervPrint (out, os.str (), comm);
  }

  // if (myRank == 0) {
  //   TEST_EQUALITY( curLineNum, expectedNumLinesRead );
  //   TEST_EQUALITY( totalNumEntriesRead, expectedTotalNumEntriesRead );
  // }
  if (myRank == 0) {
    TEST_EQUALITY( myMatrixEntries.size (), expectedNumEnt );
    if (myMatrixEntries.size () == expectedNumEnt) {
      auto iter = myMatrixEntries.begin ();
      for (std::size_t k = 0;
           k < myMatrixEntries.size () && iter != myMatrixEntries.end ();
           ++k, ++iter) {
        TEST_EQUALITY( iter->first.row, expectedRowInds[k] );
        TEST_EQUALITY( iter->first.col, expectedColInds[k] );
        TEST_EQUALITY( iter->second, expectedVals[k] );
      }
    }
  }
  else if (myRank == 1) {
    TEST_EQUALITY( myMatrixEntries.size (), expectedNumEnt );
    if (myMatrixEntries.size () == expectedNumEnt) {
      auto iter = myMatrixEntries.begin ();
      for (std::size_t k = 0;
           k < myMatrixEntries.size () && iter != myMatrixEntries.end ();
           ++k, ++iter) {
        TEST_EQUALITY( iter->first.row, expectedRowInds[k] );
        TEST_EQUALITY( iter->first.col, expectedColInds[k] );
        TEST_EQUALITY( iter->second, expectedVals[k] );
      }
    }
  }

  lclSuccess = success ? 1 : 0;
  gblSuccess = 0; // output argument
  reduceAll<int, int> (comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY( gblSuccess, 1 );
}


void
testComplex (bool& success,
             Teuchos::FancyOStream& out,
             const char testCase[],
             const GO expectedRowInds[], // on the calling process
             const GO expectedColInds[], // on the calling process
             const std::complex<double> expectedVals[], // on the calling process
             const std::size_t expectedNumEnt, // on the calling process
             const Comm<int>& comm,
             const bool debug = false,
             std::ostream* errStrm = NULL)
{
  using Tpetra::Details::readAndDealOutTriples;
  using Teuchos::outArg;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  typedef std::complex<double> SC;
  // Type representing the local entries of a sparse matrix.  The set
  // sorts its entries lexicographically, first by row index, then by
  // column index.
  typedef std::map<CooGraphEntry<GO>, SC,
                   CompareCooGraphEntries<GO> > matrix_type;
  const int myRank = comm.getRank ();
  int lclSuccess = 1;
  int gblSuccess = 0; // output argument

  out << "Test readAndDealOutTriples for complex-valued matrix entries" << endl;
  Teuchos::OSTab tab1 (out);

  // The entries of the sparse matrix that live on the calling process.
  matrix_type myMatrixEntries;
  // Closure that inserts a sparse matrix entry into myMatrixEntries
  // if it does not yet exist, else adds the input value to that of
  // the existing entry.
  std::function<int (const GO, const GO, const SC&)> insertEntry =
    [&myMatrixEntries] (const GO rowInd, const GO colInd, const SC& val) {
    auto iter = myMatrixEntries.find ({rowInd, colInd});
    if (iter == myMatrixEntries.end ()) {
      myMatrixEntries.insert ({{rowInd, colInd}, val});
    }
    else {
      iter->second += val;
    }
    return 0;
  };

  out << "Open the \"file\" on Process 0" << endl;
  std::istringstream inputStream (testCase);

  std::size_t curLineNum = 0;
  std::size_t totalNumEntriesRead = 0;

  // Keep this small, so that we can observe complete function
  // behavior.  If not less than the number of input lines, then
  // Process 0 will simply read all the entries and handle them
  // itself.  The tests below are hard-coded to assume this is 2, so
  // please change the tests if you want to change this number.
  const std::size_t maxNumEntriesPerMessage = 2;

  // One line of the input "file" is a comment line.  Note that these
  // values depend on the input "file."
  //const std::size_t expectedNumLinesRead = 8;
  //const std::size_t expectedTotalNumEntriesRead = 7;

  const bool tolerant = false;
  const int errCode =
    readAndDealOutTriples (inputStream, curLineNum, totalNumEntriesRead,
                           insertEntry, maxNumEntriesPerMessage, comm,
                           tolerant, errStrm, debug);
  TEST_EQUALITY( errCode, 0 );

  lclSuccess = success ? 1 : 0;
  gblSuccess = 0; // output argument
  reduceAll<int, int> (comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY( gblSuccess, 1 );
  if (gblSuccess != 1) {
    out << "Some process reported a nonzero error code!" << endl;
    return;
  }

  {
    std::ostringstream os;
    os << "Proc " << myRank << ": Got matrix entries: [";
    std::size_t k = 0;
    for (auto iter = myMatrixEntries.begin (); iter != myMatrixEntries.end (); ++iter, ++k) {
      os << "{row=" << iter->first.row << ", col=" << iter->first.col << ", val=" << iter->second << "}";
      if (k + 1 < myMatrixEntries.size ()) {
        os << ", ";
      }
    }
    os << "]" << endl;
    Tpetra::Details::gathervPrint (out, os.str (), comm);
  }

  // if (myRank == 0) {
  //   TEST_EQUALITY( curLineNum, expectedNumLinesRead );
  //   TEST_EQUALITY( totalNumEntriesRead, expectedTotalNumEntriesRead );
  // }
  if (myRank == 0) {
    TEST_EQUALITY( myMatrixEntries.size (), expectedNumEnt );
    if (myMatrixEntries.size () == expectedNumEnt) {
      auto iter = myMatrixEntries.begin ();
      for (std::size_t k = 0;
           k < myMatrixEntries.size () && iter != myMatrixEntries.end ();
           ++k, ++iter) {
        TEST_EQUALITY( iter->first.row, expectedRowInds[k] );
        TEST_EQUALITY( iter->first.col, expectedColInds[k] );
        TEST_EQUALITY( iter->second, expectedVals[k] );
      }
    }
  }
  else if (myRank == 1) {
    TEST_EQUALITY( myMatrixEntries.size (), expectedNumEnt );
    if (myMatrixEntries.size () == expectedNumEnt) {
      auto iter = myMatrixEntries.begin ();
      for (std::size_t k = 0;
           k < myMatrixEntries.size () && iter != myMatrixEntries.end ();
           ++k, ++iter) {
        TEST_EQUALITY( iter->first.row, expectedRowInds[k] );
        TEST_EQUALITY( iter->first.col, expectedColInds[k] );
        TEST_EQUALITY( iter->second, expectedVals[k] );
      }
    }
  }

  lclSuccess = success ? 1 : 0;
  gblSuccess = 0; // output argument
  reduceAll<int, int> (comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY( gblSuccess, 1 );
}


TEUCHOS_UNIT_TEST( ReadTriples, Real1 )
{
  auto comm = Tpetra::TestingUtilities::getDefaultComm ();
  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();

  TEST_EQUALITY( numProcs, 2 );
  if (numProcs != 2) {
    out << "This test requires exactly 2 MPI processes";
  }

#ifdef HAVE_TPETRACORE_MPI
  // Set the MPI error handler so that errors return, instead of
  // immediately causing MPI_Abort.  This will help us catch any bugs
  // with how we use, e.g., MPI_Pack and MPI_Unpack.
  {
    using Teuchos::MpiComm;
    using Teuchos::rcp_const_cast;
    using Teuchos::rcp_dynamic_cast;

    constexpr bool throwOnFail = true;
    auto mpiComm = rcp_dynamic_cast<const MpiComm<int> > (comm, throwOnFail);
    // We have to cast away const to call setErrorHandler.
    auto mpiCommNonConst = rcp_const_cast<MpiComm<int> > (mpiComm);
    auto errHandler =
      rcp (new Teuchos::OpaqueWrapper<MPI_Errhandler> (MPI_ERRORS_RETURN));
    mpiCommNonConst->setErrorHandler (errHandler);
  }
#endif // HAVE_TPETRACORE_MPI

  const bool debug = true;
  std::ostream* errStrm = debug ? &std::cerr : NULL;

  if (myRank == 0) {
    const GO expectedRowInds[] = {1, 3, 9, 11};
    const GO expectedColInds[] = {2, 4, 10, 12};
    const double expectedVals[] = {12.0, 34.0, 100.0, 122.0};
    const std::size_t expectedNumEnt = 4;
    testReal (success, out, inputFileReal,
              expectedRowInds, expectedColInds,
              expectedVals, expectedNumEnt,
              *comm, debug, errStrm);
  }
  else if (myRank == 1) {
    const GO expectedRowInds[] = {5, 7, 13};
    const GO expectedColInds[] = {6, 8, 14};
    const double expectedVals[] = {56.0, 78.0, 144.0};
    const std::size_t expectedNumEnt = 3;
    testReal (success, out, inputFileReal,
              expectedRowInds, expectedColInds,
              expectedVals, expectedNumEnt,
              *comm, debug, errStrm);
  }
  else {
    const GO* const expectedRowInds = NULL;
    const GO* const expectedColInds = NULL;
    const double* const expectedVals = NULL;
    const std::size_t expectedNumEnt = 0;
    testReal (success, out, inputFileReal,
              expectedRowInds, expectedColInds,
              expectedVals, expectedNumEnt,
              *comm, debug, errStrm);
  }
}


TEUCHOS_UNIT_TEST( ReadTriples, Real2 )
{
  auto comm = Tpetra::TestingUtilities::getDefaultComm ();
  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();

  TEST_EQUALITY( numProcs, 2 );
  if (numProcs != 2) {
    out << "This test requires exactly 2 MPI processes";
  }

#ifdef HAVE_TPETRACORE_MPI
  // Set the MPI error handler so that errors return, instead of
  // immediately causing MPI_Abort.  This will help us catch any bugs
  // with how we use, e.g., MPI_Pack and MPI_Unpack.
  {
    using Teuchos::MpiComm;
    using Teuchos::rcp_const_cast;
    using Teuchos::rcp_dynamic_cast;

    constexpr bool throwOnFail = true;
    auto mpiComm = rcp_dynamic_cast<const MpiComm<int> > (comm, throwOnFail);
    // We have to cast away const to call setErrorHandler.
    auto mpiCommNonConst = rcp_const_cast<MpiComm<int> > (mpiComm);
    auto errHandler =
      rcp (new Teuchos::OpaqueWrapper<MPI_Errhandler> (MPI_ERRORS_RETURN));
    mpiCommNonConst->setErrorHandler (errHandler);
  }
#endif // HAVE_TPETRACORE_MPI

  const bool debug = true;
  std::ostream* errStrm = debug ? &std::cerr : NULL;

  if (myRank == 0) {
    const GO expectedRowInds[] = {1, 3, 9, 11};
    const GO expectedColInds[] = {2, 4, 10, 12};
    const double expectedVals[] = {12.0, 34.0, 100.0, 122.0};
    const std::size_t expectedNumEnt = 4;
    testReal (success, out, inputFileReal2,
              expectedRowInds, expectedColInds,
              expectedVals, expectedNumEnt,
              *comm, debug, errStrm);
  }
  else if (myRank == 1) {
    const GO expectedRowInds[] = {5, 7};
    const GO expectedColInds[] = {6, 8};
    const double expectedVals[] = {56.0, 78.0};
    const std::size_t expectedNumEnt = 2;
    testReal (success, out, inputFileReal2,
              expectedRowInds, expectedColInds,
              expectedVals, expectedNumEnt,
              *comm, debug, errStrm);
  }
  else {
    const GO* const expectedRowInds = NULL;
    const GO* const expectedColInds = NULL;
    const double* const expectedVals = NULL;
    const std::size_t expectedNumEnt = 0;
    testReal (success, out, inputFileReal2,
              expectedRowInds, expectedColInds,
              expectedVals, expectedNumEnt,
              *comm, debug, errStrm);
  }
}


TEUCHOS_UNIT_TEST( ReadTriples, Real3 )
{
  auto comm = Tpetra::TestingUtilities::getDefaultComm ();
  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();

  TEST_EQUALITY( numProcs, 2 );
  if (numProcs != 2) {
    out << "This test requires exactly 2 MPI processes";
  }

#ifdef HAVE_TPETRACORE_MPI
  // Set the MPI error handler so that errors return, instead of
  // immediately causing MPI_Abort.  This will help us catch any bugs
  // with how we use, e.g., MPI_Pack and MPI_Unpack.
  {
    using Teuchos::MpiComm;
    using Teuchos::rcp_const_cast;
    using Teuchos::rcp_dynamic_cast;

    constexpr bool throwOnFail = true;
    auto mpiComm = rcp_dynamic_cast<const MpiComm<int> > (comm, throwOnFail);
    // We have to cast away const to call setErrorHandler.
    auto mpiCommNonConst = rcp_const_cast<MpiComm<int> > (mpiComm);
    auto errHandler =
      rcp (new Teuchos::OpaqueWrapper<MPI_Errhandler> (MPI_ERRORS_RETURN));
    mpiCommNonConst->setErrorHandler (errHandler);
  }
#endif // HAVE_TPETRACORE_MPI

  const bool debug = true;
  std::ostream* errStrm = debug ? &std::cerr : NULL;

  if (myRank == 0) {
    const GO expectedRowInds[] = {1, 3};
    const GO expectedColInds[] = {2, 4};
    const double expectedVals[] = {12.0, 34.0};
    const std::size_t expectedNumEnt = 2;
    testReal (success, out, inputFileReal3,
              expectedRowInds, expectedColInds,
              expectedVals, expectedNumEnt,
              *comm, debug, errStrm);
  }
  else if (myRank == 1) {
    const GO* const expectedRowInds = NULL;
    const GO* const expectedColInds = NULL;
    const double* const expectedVals = NULL;
    const std::size_t expectedNumEnt = 0;
    testReal (success, out, inputFileReal3,
              expectedRowInds, expectedColInds,
              expectedVals, expectedNumEnt,
              *comm, debug, errStrm);
  }
  else {
    const GO* const expectedRowInds = NULL;
    const GO* const expectedColInds = NULL;
    const double* const expectedVals = NULL;
    const std::size_t expectedNumEnt = 0;
    testReal (success, out, inputFileReal3,
              expectedRowInds, expectedColInds,
              expectedVals, expectedNumEnt,
              *comm, debug, errStrm);
  }
}

TEUCHOS_UNIT_TEST( ReadTriples, Complex )
{
  auto comm = Tpetra::TestingUtilities::getDefaultComm ();
  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();

  TEST_EQUALITY( numProcs, 2 );
  if (numProcs != 2) {
    out << "This test requires exactly 2 MPI processes";
  }

#ifdef HAVE_TPETRACORE_MPI
  // Set the MPI error handler so that errors return, instead of
  // immediately causing MPI_Abort.  This will help us catch any bugs
  // with how we use, e.g., MPI_Pack and MPI_Unpack.
  {
    using Teuchos::MpiComm;
    using Teuchos::rcp_const_cast;
    using Teuchos::rcp_dynamic_cast;

    constexpr bool throwOnFail = true;
    auto mpiComm = rcp_dynamic_cast<const MpiComm<int> > (comm, throwOnFail);
    // We have to cast away const to call setErrorHandler.
    auto mpiCommNonConst = rcp_const_cast<MpiComm<int> > (mpiComm);
    auto errHandler =
      rcp (new Teuchos::OpaqueWrapper<MPI_Errhandler> (MPI_ERRORS_RETURN));
    mpiCommNonConst->setErrorHandler (errHandler);
  }
#endif // HAVE_TPETRACORE_MPI

  const bool debug = true;
  std::ostream* errStrm = debug ? &std::cerr : NULL;

  if (myRank == 0) {
    const GO expectedRowInds[] = {1, 3, 9, 11};
    const GO expectedColInds[] = {2, 4, 10, 12};
    const std::complex<double> expectedVals[] = {
      {12.0, -1200.0},
      {34.0, -3400.0},
      {100.0, -10000.0},
      {122.0, -12200.0}
    };
    const std::size_t expectedNumEnt = 4;
    testComplex (success, out, inputFileComplex,
                 expectedRowInds, expectedColInds,
                 expectedVals, expectedNumEnt, *comm,
                 debug, errStrm);
  }
  else if (myRank == 1) {
    const GO expectedRowInds[] = {5, 7, 13};
    const GO expectedColInds[] = {6, 8, 14};
    const std::complex<double> expectedVals[] = {
      {56.0, -5600.0},
      {78.0, -7800.0},
      {144.0, -14400.0}
    };
    const std::size_t expectedNumEnt = 3;
    testComplex (success, out, inputFileComplex,
                 expectedRowInds, expectedColInds,
                 expectedVals, expectedNumEnt, *comm,
                 debug, errStrm);
  }
  else {
    const GO* const expectedRowInds = NULL;
    const GO* const expectedColInds = NULL;
    const std::complex<double>* const expectedVals = NULL;
    const std::size_t expectedNumEnt = 0;
    testComplex (success, out, inputFileComplex,
                 expectedRowInds, expectedColInds,
                 expectedVals, expectedNumEnt, *comm,
                 debug, errStrm);
  }
}

} // namespace (anonymous)

