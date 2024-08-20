// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <MatrixMarket_Tpetra.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_oblackholestream.hpp>

const char standardFailMessage[] = "End Result: TEST FAILED";
const char standardPassMessage[] = "End Result: TEST PASSED";

// Run the test over the given communicator.
//
// Return 1 if the test succeeded ON THE CALLING PROCESS,
// else return 0 if the test failed ON THE CALLING PROCESS.
int
testImpl (Teuchos::FancyOStream& out,
          const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::endl;
  typedef Tpetra::CrsMatrix<> crs_matrix_type;
  typedef Tpetra::Map<> map_type;
  typedef Tpetra::MultiVector<> MV;
  typedef Tpetra::MatrixMarket::Writer<crs_matrix_type> writer_type;
  typedef Tpetra::MatrixMarket::Reader<crs_matrix_type> reader_type;
  typedef MV::scalar_type scalar_type;
  typedef MV::global_ordinal_type GO;
  typedef Tpetra::global_size_t GST;

  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();

  Teuchos::OSTab tab1 (out); // indent by (another) one for the current scope

  // Bug 6288 specifically requires running on multiple MPI processes.
  if (numProcs < 2) {
    out << "This test only makes sense when run on > 1 MPI processes!"
        << endl << endl << standardFailMessage << endl;
    return 0;
  }

  const GST gblNumRows = 64;
  const size_t numCols = 2;
  {
    // Make sure we aren't running on too _many_ MPI processes.  We
    // test this by ensuring that each process has at least numCols
    // rows of the MultiVector.
    const int maxNumProcs =
      static_cast<int> (gblNumRows) / static_cast<int> (numCols);
    if (numProcs > maxNumProcs) {
      out << "This test only makes sense when run on <= " << maxNumProcs
          << " MPI processes!" << endl << endl << standardFailMessage << endl;
      return 0;
    }
  }

  const GO indexBase = 0;
  RCP<const map_type> origMap = rcp (new map_type (gblNumRows, indexBase, comm));

  out << "Create MultiVector" << endl;

  RCP<MV> mv = rcp (new MV (origMap, numCols));
  for (GO i = origMap->getMinGlobalIndex ();
       i < origMap->getMaxGlobalIndex () + 1; ++i) {
    for (size_t j = 0; j < numCols; ++j) {
      mv->replaceGlobalValue (i, j, j*64 + i);
    }
  }

  out << "Write MultiVector and Map" << endl;

  std::ostringstream mvOutFile;
  writer_type::writeDense (mvOutFile, *mv);
  std::ostringstream mapOutFile;
  writer_type::writeMap (mapOutFile, * (mv->getMap ()));

  out << "Read Map and MultiVector" << endl;

  std::istringstream mapInFile (mapOutFile.str ());
  RCP<const map_type> readMap =
    reader_type::readMap (mapInFile, comm,
                          false, false);
  std::istringstream mvInFile (mvOutFile.str ());
  RCP<MV> read_mv = reader_type::readDense (mvInFile, comm, 
                                            readMap, false, false);

  out << "Write MultiVector again" << endl;

  std::ostringstream mvOutFile2;
  writer_type::writeDense (mvOutFile2, read_mv);

  out << "Test expected output against both actual output instances" << endl;

  int lclSuccess = 1;
  if (myRank == 0) {
    std::ostringstream exp;

    // FIXME (mfh 22 Jul 2015) This forces scientific notation and a
    // minimum number of digits to represent the precision correctly.
    // The test really shouldn't depend as much on the printed digits;
    // it should depend on the values not changing.
    Teuchos::SetScientific<scalar_type> sci (exp);

    exp << "%%MatrixMarket matrix array real general" << endl
        << gblNumRows << " " << numCols << endl;
    for (size_t j = 0; j < numCols; ++j) {
      for (size_t i = 0; i < origMap->getGlobalNumElements (); ++i) {
        const scalar_type val = static_cast<scalar_type> (j*64 + i);
        exp << val << endl;
      }
    }
    const std::string expStr = exp.str ();
    const bool firstEqual = (expStr == mvOutFile.str ());
    const bool secondEqual = (expStr == mvOutFile2.str ());

    if (! firstEqual) {
      out << "First comparison FAILED" << endl;
    }
    if (! secondEqual) {
      out << "Second comparison FAILED" << endl;
    }

    lclSuccess = (firstEqual && secondEqual) ? 1 : 0;
  }
  return lclSuccess;
}

// Run the test over the given communicator.
//
// Return 1 GLOBALLY if the test succeeded on ALL PROCESSES in the
// given communicator, else return 0 GLOBALLY.
int
test (const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
{
  using Teuchos::broadcast;
  using Teuchos::inOutArg;
  using Teuchos::RCP;
  using std::endl;

  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();

  // Make an output stream (for verbose output) that only prints on
  // Process 0 of the given communicator.  "FancyOStream" is an
  // std::ostream for which we can define automatic indentation.
  Teuchos::oblackholestream blackHole;
  std::ostream& outStrm = (myRank == 0) ? std::cout : blackHole;
  RCP<Teuchos::FancyOStream> fancyOutStrmPtr =
    Teuchos::getFancyOStream (Teuchos::rcpFromRef (outStrm));
  Teuchos::FancyOStream& out = *fancyOutStrmPtr;

  Teuchos::OSTab tab0 (out); // indent by one for the current scope
  out << "Testing on " << numProcs << " MPI process(es)" << endl;

  int lclSuccess = testImpl (out, comm);
  int gblSuccess = lclSuccess;
  broadcast<int, int> (*comm, 0, inOutArg (gblSuccess));
  out << (gblSuccess == 1 ? standardPassMessage : standardFailMessage) << endl;

  return gblSuccess;
}


int
main (int argc, char **argv)
{
  using Teuchos::broadcast;
  using Teuchos::Comm;
  using Teuchos::MpiComm;
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::endl;
  // Common text in exception messages below.
  const char suffix[] = ".  This suggests a bug in Teuchos::MpiComm::split.  "
    "Please report this bug to the Teuchos and Tpetra developers.";

  int curGblSuccess = 1;
  int gblSuccess = 1;

  MPI_Init (&argc, &argv);
  RCP<const Comm<int> > commWorld = rcp (new MpiComm<int> (MPI_COMM_WORLD));
  const int myRank = commWorld->getRank ();
  const int numProcs = commWorld->getSize ();

  Teuchos::oblackholestream blackHole;
  std::ostream& out = (commWorld->getRank () == 0) ? std::cout : blackHole;
  out << "Test for Bug 6288 (Tpetra::MultiVector Matrix Market I/O)" << endl;

  // Test over the full communicator.
  curGblSuccess = test (commWorld);
  gblSuccess = (curGblSuccess == 1 && gblSuccess == 1) ? 1 : 0;

  // Repeat test for size >= 2 subsets of the full communicator.  We
  // include the newNumProcs == numProcs case as a sanity check for
  // Comm::split().
  for (int newNumProcs = 2; newNumProcs <= numProcs; newNumProcs *= 2) {
    // Split commWorld into two communicators.  Color 0 gets the
    // processes that will participate in the current test.  Color 1
    // gets all the other processes (that don't participate).
    const int color = (myRank < newNumProcs) ? 0 : 1;
    const int key = 0; // we only need to use the 'color' argument here
    RCP<Comm<int> > curComm = commWorld->split (color, key);

    TEUCHOS_TEST_FOR_EXCEPTION
      (color == 0 && curComm->getSize () != newNumProcs, std::logic_error,
       "For color = " << color << ", curComm->getSize() = " <<
       curComm->getSize () << " != newNumProcs = " << newNumProcs << suffix);
    TEUCHOS_TEST_FOR_EXCEPTION
      (color == 1 && curComm->getSize () != (numProcs - newNumProcs),
       std::logic_error, "For color = " << color << ", curComm->getSize() = "
       << curComm->getSize () << " != numProcs - newNumProcs = "
       << (numProcs - newNumProcs) << suffix);

    if (color == 0) { // participating processes
      curGblSuccess = test (curComm);
      gblSuccess = (curGblSuccess == 1 && gblSuccess == 1) ? 1 : 0;
    }
  }

  // The above tests over subsets of commWorld only do a "global"
  // check over the subset, not over commWorld.  Thus, we do a final
  // check here over commWorld.
  broadcast<int, int> (*commWorld, 0, outArg (gblSuccess));
  out << (gblSuccess == 1 ? standardPassMessage : standardFailMessage) << endl;

  MPI_Finalize ();
  return (gblSuccess == 1) ? EXIT_SUCCESS : EXIT_FAILURE;
}
