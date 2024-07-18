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
#include "Kokkos_Core.hpp"

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

int getRankInCommWorld ()
{
  int myRank = 0;
  (void) MPI_Comm_rank (MPI_COMM_WORLD, &myRank);
  return myRank;
}

bool allTrueInCommWorld (const bool lclTruth)
{
  int lclTruthInt = lclTruth ? 1 : 0;
  int gblTruthInt = 0;
  MPI_Allreduce (&lclTruthInt, &gblTruthInt, 1, MPI_INT,
                 MPI_MIN, MPI_COMM_WORLD);
  return gblTruthInt != 0;
}

bool
tpetraCommIsLocallyLegit (const Teuchos::Comm<int>* wrappedTpetraComm)
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
  (void) MPI_Comm_compare (MPI_COMM_WORLD, tpetraComm, &result);
  // Tpetra reserves the right to MPI_Comm_dup on the input comm.
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
    cout << "MPI_Initialized claims MPI is initialized, "
      "before MPI_Init was called" << endl;
    return;
  }
  (void) MPI_Init (&argc, &argv);
  if (! isMpiInitialized ()) {
    success = false;
    cout << "MPI_Initialized claims MPI is not initialized, "
      "even after MPI_Init was called" << endl;
    return;
  }
  const int myRank = getRankInCommWorld ();

  if (Kokkos::is_initialized ()) {
    success = false;
    cout << "Kokkos::is_initialized() is true, "
      "even though Kokkos::initialize was not yet called." << endl;
    return;
  }

  // This test exercises the case where the "user" has called MPI_Init
  // before calling Tpetra::ScopeGuard's constructor, but has not yet
  // called Kokkos::initialize.
  {
    Tpetra::ScopeGuard tpetraScope (&argc, &argv);
    if (! isMpiInitialized ()) {
      success = false;
      cout << "MPI_Initialized claims MPI was not initialized, "
        "even after MPI_Init and Tpetra::ScopeGuard::ScopeGuard were called."
        << endl;
      return;
    }
    if (! Kokkos::is_initialized ()) {
      success = false;
      cout << "Kokkos::is_initialized() is false, "
        "even after Tpetra::ScopeGuard::ScopeGuard was called." << endl;
      return;
    }

    // MPI is initialized, so we can check whether all processes
    // report Tpetra as initialized.
    const bool tpetraIsNowInitialized =
      allTrueInCommWorld (Tpetra::isInitialized ());
    if (! tpetraIsNowInitialized) {
      success = false;
      if (myRank == 0) {
        cout << "Tpetra::isInitialized() is false on at least one process"
          ", even after Tpetra::ScopeGuard::ScopeGuard was called." << endl;
      }
      if (Kokkos::is_initialized ()) {
        Kokkos::finalize ();
      }
      (void) MPI_Finalize (); // just for completeness
      return;
    }

    auto comm = Tpetra::getDefaultComm ();
    const bool tpetraCommGloballyValid =
      allTrueInCommWorld (tpetraCommIsLocallyLegit (comm.get ()));
    if (! tpetraCommGloballyValid) {
      success = false;
      if (myRank == 0) {
        cout << "Tpetra::getDefaultComm() returns an invalid comm "
          "on at least one process." << endl;
      }
    }

    const int myTpetraRank = comm.is_null () ? 0 : comm->getRank ();
    const bool ranksSame = allTrueInCommWorld (myRank == myTpetraRank);
    if (! ranksSame) {
      success = false;
      if (myRank == 0) {
        cout << "MPI rank does not match Tpetra rank "
          "on at least one process" << endl;
      }
    }

    if (myRank == 0) {
      cout << "About to leave Tpetra scope" << endl;
    }
  }
  if (myRank == 0) {
    cout << "Left Tpetra scope" << endl;
  }
  // Tpetra::ScopeGuard's destructor was responsible for calling
  // Kokkos::finalize.
  if (Kokkos::is_initialized ()) {
    success = false;
    cout << "Kokkos::is_initialized() is true, "
      "after Tpetra::ScopeGuard::~ScopeGuard was called." << endl;
  }
  // Since the "user" is responsible for calling MPI_Finalize,
  // Tpetra::ScopeGuard's destructor should NOT have called
  // MPI_Finalize.
  if (! isMpiInitialized ()) {
    success = false;
    cout << "Tpetra::ScopeGuard::~ScopeGuard seems to have called "
      "MPI_Finalize, even though the user was responsible for initializing "
      "and finalizing MPI." << endl;
    return;
  }

  // MPI is still initialized, so we can check whether processes are
  // consistent.
  const bool tpetraGloballyFinalized =
    allTrueInCommWorld (! Tpetra::isInitialized ());
  if (! tpetraGloballyFinalized) {
    success = false;
    if (myRank == 0) {
      cout << "Tpetra::isInitialized() returns true on some process, "
        "even after Tpetra::ScopeGuard::~ScopeGuard has been called." << endl;
    }
  }

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
