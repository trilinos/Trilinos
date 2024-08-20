// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <cstdlib>
#include <iostream>
#include "Tpetra_Core.hpp"

#if ! defined(HAVE_TPETRACORE_MPI)
#  error "Building and testing this example requires MPI."
#endif // ! defined(HAVE_TPETRACORE_MPI)
#include "mpi.h"
#include "Tpetra_Details_extractMpiCommFromTeuchos.hpp"


bool isMpiInitialized ()
{
  int mpiInitializedInt = 0;
  (void) MPI_Initialized (&mpiInitializedInt);
  return mpiInitializedInt != 0;
}

int getRankInComm (MPI_Comm comm)
{
  int myRank = 0;
  (void) MPI_Comm_rank (comm, &myRank);
  return myRank;
}

bool allTrueInComm (const bool lclTruth, MPI_Comm comm)
{
  int lclTruthInt = lclTruth ? 1 : 0;
  int gblTruthInt = 0;
  (void) MPI_Allreduce (&lclTruthInt, &gblTruthInt, 1,
                        MPI_INT, MPI_MIN, comm);
  return gblTruthInt != 0;
}

bool
tpetraCommIsLocallyLegit (const Teuchos::Comm<int>* wrappedTpetraComm,
                          MPI_Comm expectedMpiComm)
{
  if (wrappedTpetraComm == nullptr) {
    return false;
  }
  MPI_Comm tpetraComm;
  try {
    using Tpetra::Details::extractMpiCommFromTeuchos;
    tpetraComm = extractMpiCommFromTeuchos (*wrappedTpetraComm);
  }
  catch (...) {
    return false;
  }
  if (tpetraComm == MPI_COMM_NULL) {
    return false;
  }
  int result = MPI_UNEQUAL;
  (void) MPI_Comm_compare (expectedMpiComm, tpetraComm, &result);
  // Tpetra reserves the right to MPI_Comm_dup the input comm.
  return result == MPI_IDENT || result == MPI_CONGRUENT;
}


// NOTE TO TEST AUTHORS: The code that calls this function captures
// std::cerr, so don't write to std::cerr on purpose in this function.
void testMain (bool& success, int argc, char* argv[])
{
  using std::cout;
  using std::endl;

  if (isMpiInitialized ()) {
    success = false;
    cout << "TEST FAILED: MPI_Initialized claims MPI is "
      "initialized, before MPI_Init was called." << endl;
    return;
  }
  (void) MPI_Init (&argc, &argv);
  if (! isMpiInitialized ()) {
    success = false;
    cout << "TEST FAILED: MPI_Initialized claims MPI is not "
      "initialized, even after MPI_Init was called." << endl;
    return;
  }

  // Split off Process 0 into a comm by itself.  Only create a
  // Tpetra::ScopeGuard instance on the remaining processes.  The
  // point is to test that there are no hidden dependencies on
  // MPI_COMM_WORLD, since those often manifest as asking Process 0 in
  // MPI_COMM_WORLD to do something.
  const int color = (getRankInComm (MPI_COMM_WORLD) == 0) ? 0 : 1;
  MPI_Comm splitComm = MPI_COMM_NULL;
  {
    const int key = 0;
    const int errCode =
      MPI_Comm_split (MPI_COMM_WORLD, color, key, &splitComm);
    if (errCode != MPI_SUCCESS) {
      success = false;
      cout << "TEST FAILED: MPI_Comm_split failed!" << endl;
      (void) MPI_Abort (MPI_COMM_WORLD, EXIT_FAILURE);
    }
  }

  if (color == 0) { // this subcomm doesn't participate
    MPI_Comm_free (&splitComm);
    MPI_Finalize ();
    return;
  }
  // Invariant: color == 1, and we're on the corresponding splitComm.
  const int mySplitRank = getRankInComm (splitComm);

  // Special case of Tpetra::ScopeGuard, for when users want Tpetra to
  // use a custom MPI_Comm as Tpetra's default communicator.
  {
    Tpetra::ScopeGuard tpetraScope (&argc, &argv, splitComm);
    if (! isMpiInitialized ()) {
      success = false;
      cout << "TEST FAILED: MPI_Initialized claims MPI was not initialized, "
        "even after MPI_Init and Tpetra::ScopeGuard::ScopeGuard were called."
        << endl;
      return;
    }

    // MPI is initialized, so we can check whether all processes
    // report Tpetra as initialized.
    const bool tpetraIsNowInitialized =
      allTrueInComm (Tpetra::isInitialized (), splitComm);
    if (! tpetraIsNowInitialized) {
      success = false;
      if (mySplitRank == 0) {
        cout << "TEST FAILED: Tpetra::isInitialized() is false on at least "
          "one process, even after Tpetra::ScopeGuard::ScopeGuard was called."
          << endl;
      }
      (void) MPI_Abort (MPI_COMM_WORLD, EXIT_FAILURE);
    }

    auto comm = Tpetra::getDefaultComm ();
    const bool tpetraCommLocallyValid =
      tpetraCommIsLocallyLegit (comm.get (), splitComm);
    const bool tpetraCommGloballyValid =
      allTrueInComm (tpetraCommLocallyValid, splitComm);
    if (! tpetraCommGloballyValid) {
      success = false;
      if (mySplitRank == 0) {
        cout << "TEST FAILED: Tpetra::getDefaultComm() returns "
          "an invalid comm on at least one process." << endl;
      }
    }

    {
      // The above check already tests whether comm is null, so we
      // just need to make sure that error case won't dereference
      // nullptr.
      const int myTpetraRank = comm.is_null () ? 0 : comm->getRank ();
      const bool ranksSame =
        allTrueInComm (mySplitRank == myTpetraRank, splitComm);
      if (! ranksSame) {
        success = false;
        if (mySplitRank == 0) {
          cout << "TEST FAILED: MPI rank does not match Tpetra rank "
            "on at least one process." << endl;
        }
      }
    }

    if (mySplitRank == 0) {
      cout << "About to leave Tpetra scope" << endl;
    }
  }
  if (mySplitRank == 0) {
    cout << "Left Tpetra scope" << endl;
  }
  // Since the "user" is responsible for calling MPI_Finalize,
  // Tpetra::ScopeGuard's destructor should NOT have called
  // MPI_Finalize.
  if (! isMpiInitialized ()) {
    success = false;
    cout << "TEST FAILED: Tpetra::ScopeGuard::~ScopeGuard seems to have "
      "called MPI_Finalize, even though the user was supposed to be "
      "responsible for initializing and finalizing MPI." << endl;
    return;
  }

  // MPI is still initialized, so we can check whether processes are
  // consistent.
  const bool tpetraGloballyFinalized =
    allTrueInComm (! Tpetra::isInitialized (), splitComm);
  if (! tpetraGloballyFinalized) {
    success = false;
    if (mySplitRank == 0) {
      cout << "TEST FAILED: Tpetra::isInitialized() returns true on some "
        "process, even after Tpetra::ScopeGuard::~ScopeGuard has been called."
        << endl;
    }
  }

  (void) MPI_Comm_free (&splitComm);
  (void) MPI_Finalize ();
}

int main (int argc, char* argv[])
{
  using std::cout;
  using std::endl;

  bool success = true;
  testMain (success, argc, argv);

  cout << "End Result: TEST " << (success ? "PASSED" : "FAILED") << endl;
  return EXIT_SUCCESS;
}
