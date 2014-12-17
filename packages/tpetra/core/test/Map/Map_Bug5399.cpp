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

// Macros to ensure that if CUDA and KokkosCompat are enabled, then
// only the .cu (CUDA) version of this file will be compiled.
#include <Tpetra_config.h>

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Tpetra_Map.hpp>

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
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
  typedef long long GO;
  if (sizeof (long long) <= 4) {
    out << "sizeof (long long) = " << sizeof (long long) << " <= 4.  "
        << "This test only makes sense if sizeof (long long) >= 8." << endl;
    return;
  }
#else // NOT HAVE_TEUCHOS_LONG_LONG_INT
  typedef long GO;
  if (sizeof (long) <= 4) {
    out << "sizeof (long) = " << sizeof (long) << " <= 4.  "
        << "This test only makes sense if sizeof (long) >= 8." << endl;
    return;
  }
#endif // HAVE_TEUCHOS_LONG_LONG_INT
  typedef Tpetra::Map<int, GO> map_type;

  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm ();
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

