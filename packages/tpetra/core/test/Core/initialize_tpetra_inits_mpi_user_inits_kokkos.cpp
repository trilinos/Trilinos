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
#include <stdexcept>
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

bool isMpiFinalized ()
{
  int mpiFinalizedInt = 0;
  (void) MPI_Finalized (&mpiFinalizedInt);
  return mpiFinalizedInt != 0;
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

  // In this example, Tpetra::initialize is responsible for calling
  // MPI_Init and MPI_Finalize, but the "user" is responsible for
  // calling Kokkos::initialize and Kokkos::finalize.
  if (isMpiInitialized ()) {
    success = false;
    cout << "MPI_Initialized claims MPI was initialized, "
      "before Tpetra::initialize was called." << endl;
    return;
  }

  if (Kokkos::is_initialized ()) {
    success = false;
    cout << "Kokkos::is_initialized() is true, "
      "before Kokkos::initialize was called." << endl;
    return;
  }

  { // GlobalMPISession is in local scope so we can check MPI_Finalize

  Teuchos::GlobalMPISession mpiSession(&argc, &argv); // before Kokkos::initialize
  Kokkos::initialize (argc, argv);
  if (! Kokkos::is_initialized ()) {
    success = false;
    cout << "Kokkos::is_initialized() is false, "
      "after Kokkos::initialize was called." << endl;
    return;
  }

  cout << "Calling Tpetra::initialize" << endl;
  try {
    Tpetra::initialize (&argc, &argv);
  }
  catch (std::exception& e) {
    success = false;
    cout << "Tpetra::initialize threw an exception: "
	 << e.what () << endl;
  }
  catch (...) {
    success = false;
    cout << "Tpetra::initialize threw an exception "
      "not a subclass of std::exception." << endl;
  }
  cout << "Done with Tpetra::initialize" << endl;
  
  if (! isMpiInitialized ()) {
    success = false;
    cout << "MPI_Initialized claims MPI is not initialized, "
      "even after MPI_Init and Tpetra::initialize were called." << endl;
    Kokkos::finalize (); // for completeness
    return;
  }
  // MPI claims to have been initialized, so after this point and
  // before Tpetra::finalize, it should be safe to use MPI functions.
  
  const int myRank = getRankInCommWorld ();

  // Do all processes report Tpetra as initialized?
  const bool tpetraInitialized = allTrueInCommWorld (Tpetra::isInitialized ());
  if (! tpetraInitialized) {
    success = false;
    if (myRank == 0) {
      cout << "Tpetra::isInitialized() is false on at least one process, "
	"even after Tpetra::initialize() has been called." << endl;
    }
  }
  
  if (! Kokkos::is_initialized ()) {
    success = false;
    cout << "Kokkos::is_initialized() is false, even after "
      "Kokkos::initialize and Tpetra::initialize were called." << endl;
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

  // Means we should run on at least two MPI processes.
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
    cout << "About to call Tpetra::finalize" << endl;
  }
  Tpetra::finalize ();
  if (myRank == 0) {
    cout << "Called Tpetra::finalize" << endl;
  }

  // Tpetra::finalize is not responsible for calling
  // MPI_Finalize because GlobalMPISession was called first.
  if (isMpiFinalized ()) {
    success = false;
    cout << "Tpetra::finalize called MPI_Finalize "
         << "but should not have since GlobalMPISession controlled it."
         << endl;
  }

  // Kokkos is like Tpetra; Kokkos::is_initialized() means "was
  // initialized and was not finalized."  That differs from MPI, where
  // MPI_Initialized only refers to MPI_Init and MPI_Finalized only
  // refers to MPI_Finalize.
  if (! Kokkos::is_initialized ()) {
    success = false;
    cout << "Tpetra::finalize was not supposed to call "
      "Kokkos::finalize, but apparently did." << endl;
    return;
  }

  // MPI is no longer initialized, so we can't all-reduce on this.
  if (Tpetra::isInitialized ()) {
    success = false;
    if (myRank == 0) {
      cout << "Tpetra::isInitialized() returns true, "
	"even after Tpetra::finalize() has been called" << endl;
    }
  }

  Kokkos::finalize ();
  if (Kokkos::is_initialized ()) {
    success = false;
    cout << "Kokkos::is_initialized() is false "
      "even after calling Kokkos::finalize." << endl;
  }

  } // end of GlobalMPISession scope

  // Since GlobalMPISession's destructor was responsible for calling
  // MPI_Finalize, the destructor MUST have called MPI_Finalize.
  if (! isMpiFinalized ()) {
    success = false;
    cout << "GlobalMPISession::~GlobalMPISession did not call MPI_Finalize."
         << endl;
  }
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
