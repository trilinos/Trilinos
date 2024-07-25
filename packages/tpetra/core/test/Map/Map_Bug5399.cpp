// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Teuchos_CommHelpers.hpp"

// Unit tests like to be in an anonymous namespace.
namespace {

//
// This test must be run with 4 MPI processes.
// It ensures that Bug 5399 won't come back.
// Thanks to Lee Ann Riesen for the test case.
//
TEUCHOS_UNIT_TEST( Map, Bug5399 )
{
  using Teuchos::Array;
  using Teuchos::Comm;
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  // We use cerr rather than out for reporting exceptions on different
  // MPI processes, since the 'out' output stream provided by the
  // Teuchos unit test framework only prints on MPI process 0.
  using std::cerr;
  using std::endl;
#ifdef HAVE_TPETRA_INT_LONG_LONG
  // C++11 guarantees that sizeof(long long) >= 8.
  using GO = long long;
#else // NOT HAVE_TPETRA_INT_LONG_LONG
  using GO = long;
  // long is 32 bits on some platforms, including Windows.
  if (sizeof (long) <= 4) {
    out << "sizeof (long) = " << sizeof (long) << " <= 4.  "
        << "This test only makes sense if sizeof (long) >= 8." << endl;
    return;
  }
#endif // HAVE_TPETRA_INT_LONG_LONG
  using LO = Tpetra::Map<>::local_ordinal_type;
  using map_type = Tpetra::Map<LO, GO>;

  RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();

  TEUCHOS_TEST_FOR_EXCEPTION(
    numProcs != 4, std::logic_error, "This Bug 5399 test only works if the "
    "number of MPI processes is exactly 4.  If you are running this test "
    "manually (by calling mpirun or mpiexec with the test executable), please "
    "use exactly 4 MPI processes.  If you see this message when running using "
    "the Trilinos test framework (or via CTest), then there might be an error "
    "in the CMakeLists.txt file in this directory.  Please change that file "
    "so that it runs this test only in an MPI build, and only with exactly 4 "
    "processes.");

  const size_t nrows = (myRank == 0) ? 7 : 6;
  const size_t ngrows = 25;

  Array<GO> smallConsecutive (nrows);
  Array<GO> smallNonConsecutive (nrows);
  Array<GO> largeConsecutive (nrows);
  Array<GO> largeNonConsecutive (nrows);

  const GO offset = 0x70f000000000;

  if (myRank == 0) {
    smallConsecutive[0]=2;
    smallConsecutive[1]=3;
    smallConsecutive[2]=4;
    smallConsecutive[3]=6;
    smallConsecutive[4]=7;
    smallConsecutive[5]=8;
    smallConsecutive[6]=9;
  }

  if (myRank == 1) {
    smallConsecutive[0]=10;
    smallConsecutive[1]=11;
    smallConsecutive[2]=12;
    smallConsecutive[3]=13;
    smallConsecutive[4]=14;
    smallConsecutive[5]=15;
  }

  if (myRank == 2) {
    smallConsecutive[0]=16;
    smallConsecutive[1]=17;
    smallConsecutive[2]=18;
    smallConsecutive[3]=19;
    smallConsecutive[4]=20;
    smallConsecutive[5]=21;
  }

  if (myRank == 3) {
    smallConsecutive[0]=22;
    smallConsecutive[1]=23;
    smallConsecutive[2]=24;
    smallConsecutive[3]=25;
    smallConsecutive[4]=26;
    smallConsecutive[5]=27;
  }

  for (size_t i = 0; i < nrows; ++i) {
    largeConsecutive[i] = offset + smallConsecutive[i];
  }

  //
  // First test with small consecutive IDs
  //
  {
    int lclSuccess = 1;
    std::ostringstream os; // for exception test output on each process
    try {
      RCP<map_type> smallConsecutiveMap =
        rcp (new map_type (ngrows, smallConsecutive, 0, comm));
    }
    catch (std::exception &e) {
      lclSuccess = 0;
      os << "Proc " << myRank << ": First test with small consecutive IDs "
         << "failed: " << e.what () << endl;
    }
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    if (! gblSuccess) {
      if (myRank == 0) {
        cerr << "First test with small consecutive IDs failed on at least one "
             << "process.  Here are the exception messages (if any) on all "
             << "processes:" << endl;
      }
      for (int p = 0; p < numProcs; ++p) {
        if (myRank == p) {
          cerr << os.str () << endl;
        }
        comm->barrier (); // give output time to finish
        comm->barrier ();
        comm->barrier ();
      }
    }
  }

  //
  // First test with large consecutive IDs
  //
  {
    int lclSuccess = 1;
    std::ostringstream os; // for exception test output on each process
    try {
      RCP<map_type> largeConsecutiveMap =
        rcp (new map_type (ngrows, largeConsecutive, 0, comm));
    }
    catch (std::exception &e) {
      lclSuccess = 0;
      os << "Proc " << myRank << ": First test with large consecutive IDs "
         << "failed: " << e.what () << endl;
    }
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    if (! gblSuccess) {
      if (myRank == 0) {
        cerr << "First test with large consecutive IDs failed on at least one "
             << "process.  Here are the exception messages (if any) on all "
             << "processes:" << endl;
      }
      for (int p = 0; p < numProcs; ++p) {
        if (myRank == p) {
          cerr << os.str () << endl;
        }
        comm->barrier (); // give output time to finish
        comm->barrier ();
        comm->barrier ();
      }
    }
  }

  smallNonConsecutive[0] = smallConsecutive[0] + 28;
  largeNonConsecutive[0] = largeConsecutive[0] + 28;

  //
  // Second test with small consecutive IDs
  //
  {
    int lclSuccess = 1;
    std::ostringstream os; // for exception test output on each process
    try {
      RCP<map_type> smallConsecutiveMap =
        rcp (new map_type (ngrows, smallConsecutive, 0, comm));
    }
    catch (std::exception &e) {
      lclSuccess = 0;
      os << "Proc " << myRank << ": Second test with small consecutive IDs "
         << "failed: " << e.what () << endl;
    }
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    if (! gblSuccess) {
      if (myRank == 0) {
        cerr << "Second test with small consecutive IDs failed on at least "
             << "one process.  Here are the exception messages (if any) on "
             << "all processes:" << endl;
      }
      for (int p = 0; p < numProcs; ++p) {
        if (myRank == p) {
          cerr << os.str () << endl;
        }
        comm->barrier (); // give output time to finish
        comm->barrier ();
        comm->barrier ();
      }
    }
  }

  //
  // Second test with large consecutive IDs
  //
  {
    int lclSuccess = 1;
    std::ostringstream os; // for exception test output on each process
    try {
      RCP<map_type> largeConsecutiveMap =
        rcp (new map_type (ngrows, largeConsecutive, 0, comm));
    }
    catch (std::exception &e) {
      lclSuccess = 0;
      os << "Proc " << myRank << ": Second test with large consecutive IDs "
         << "failed: " << e.what () << endl;
    }
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );

    if (! gblSuccess) {
      if (myRank == 0) {
        cerr << "Second test with large consecutive IDs failed on at least "
             << "one process.  Here are the exception messages (if any) on "
             << "all processes:" << endl;
      }
      for (int p = 0; p < numProcs; ++p) {
        if (myRank == p) {
          cerr << os.str () << endl;
        }
        comm->barrier (); // give output time to finish
        comm->barrier ();
        comm->barrier ();
      }
    }
  }
}

} // namespace (anonymous)

