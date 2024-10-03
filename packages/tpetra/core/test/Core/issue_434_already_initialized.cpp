// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Test the case where Tpetra::initialize is called, when both Kokkos
// and MPI have been initialized.  Tpetra::initialize must, in that
// case, NOT attempt to initialize either Kokkos or MPI.  Furthermore,
// Tpetra::finalize NOT attempt to finalize either Kokkos or MPI.

#include "Tpetra_Core.hpp"
#include "Kokkos_Core.hpp"
#ifdef HAVE_TPETRACORE_MPI
#  include "mpi.h" // need direct MPI calls
#endif // HAVE_TPETRACORE_MPI
#include <cstdlib> // EXIT_SUCCESS, EXIT_FAILURE
#include <iostream>

  // Is Kokkos initialized?
  bool kokkosInitialized ()
  {
    // Kokkos lacks a global is_initialized function.  However,
    // Kokkos::initialize always initializes the default execution
    // space, so we can check that.
    return Kokkos::is_initialized ();
  }

#ifdef HAVE_TPETRACORE_MPI
  // Has MPI_Init been called?
  bool mpiInitialized ()
  {
    int mpiIsInitialized = 0;
    int mpiErr = MPI_Initialized (&mpiIsInitialized);

    if (mpiErr != MPI_SUCCESS) {
      std::cerr << "*** MPI_Initialized returned err != MPI_SUCCESS"
                << std::endl;
      // If MPI_Initialized fails, MPI is totally messed up.
      return false;
    }
    else {
      return (mpiIsInitialized != 0);
    }
  }

  // Has MPI_Finalize been called?
  bool mpiFinalized ()
  {
    int mpiIsFinalized = 0;
    int mpiErr = MPI_Finalized (&mpiIsFinalized);

    if (mpiErr != MPI_SUCCESS) {
      std::cerr << "*** MPI_Finalized returned err != MPI_SUCCESS"
                << std::endl;
      // If MPI_Finalized fails, MPI is totally messed up.
      return false;
    }
    else {
      return (mpiIsFinalized != 0);
    }
  }
#endif // HAVE_TPETRACORE_MPI

int main (int argc, char* argv[])
{
  using std::endl;
  bool success = true;
  //std::ostream& out = std::cout;
  std::ostream& err = std::cerr;

  // Initialize MPI (if applicable).
#ifdef HAVE_TPETRACORE_MPI
  {
    const int mpiErr = MPI_Init (&argc, &argv);
    if (mpiErr != MPI_SUCCESS) {
      err << "*** MPI_Init returned err != MPI_SUCCESS" << endl;
      success = false;
      // If MPI_Init fails, then the best we can do is stop the program.
      goto EndOfTest;
    }
  }
#endif // HAVE_TPETRACORE_MPI

  Kokkos::initialize (argc, argv);

  // Both Kokkos and MPI have been initialized.  Thus,
  // Tpetra::initialize must not attempt to initialize either of them.

  Tpetra::initialize (&argc, &argv);

  // MPI still needs to be initialized at this point, even though
  // Tpetra::initialize wasn't supposed to have done it.
#ifdef HAVE_TPETRACORE_MPI
  {
    int mpiIsInitialized = 0;
    int mpiErr = MPI_Initialized (&mpiIsInitialized);

    if (mpiErr != MPI_SUCCESS) {
      err << "*** MPI_Initialized returned err != MPI_SUCCESS" << endl;
      success = false;
      // If MPI_Initialized fails, MPI is totally messed up.  The best
      // we can do at this point is stop the program.
      goto EndOfTest;
    }
    else {
      if (mpiIsInitialized == 0) {
        err << "*** MPI_Initialized says MPI was not initialized" << endl;
        success = false;
        // If MPI is not initialized, it's safe just to stop the program.
        goto EndOfTest;
      }
    }
  }
#endif // HAVE_TPETRACORE_MPI

  // Kokkos still needs to be initialized at this point, even though
  // Tpetra::initialize wasn't supposed to have done it.
  if (! kokkosInitialized ()) {
    err << "*** Kokkos is not initialized" << endl;
    success = false;
  }

  {
    // Must hide this from the label, else the compiler reports an error.
    Teuchos::RCP<const Teuchos::Comm<int> > comm;
    bool threw = false;
    try {
      comm = Tpetra::getDefaultComm ();
    }
    catch (std::exception& e) {
      err << "*** Tpetra::getDefaultComm threw an std::exception: "
          << e.what () << endl;
      threw = true;
      success = false;
    }
    catch (...) {
      err << "*** Tpetra::getDefaultComm threw an exception not a subclass of "
        "std::exception" << endl;
      threw = true;
      success = false;
    }

    if (! threw && comm.is_null ()) {
      err << "*** Tpetra::getDefaultComm returned null (but did not throw an "
        "exception) after Tpetra was initialized" << endl;
      success = false;
    }
  }

  Tpetra::finalize ();

  // Tpetra::finalize is NOT supposed to finalize Kokkos or MPI,
  // because they were both initialized when Tpetra::initialize was
  // called.

#ifdef HAVE_TPETRACORE_MPI
  {
    int mpiIsFinalized = 0;
    int mpiErr = MPI_Finalized (&mpiIsFinalized);

    if (mpiErr != MPI_SUCCESS) {
      err << "*** MPI_Finalized returned err != MPI_SUCCESS" << endl;
      success = false;
      // If MPI_Finalized fails, MPI is totally messed up.  The best
      // we can do at this point is stop the program.
      goto EndOfTest;
    }
    else {
      if (mpiIsFinalized == 1) {
        err << "*** MPI_Finalized says MPI was finalized" << endl;
        success = false;
        goto EndOfTest; // MPI is finalized, so this is safe to do
      }
    }
  }
#endif // HAVE_TPETRACORE_MPI

  // Kokkos still needs to be initialized at this point, even though
  // Tpetra::initialize wasn't supposed to have done it.
  if (! kokkosInitialized ()) {
    err << "*** Kokkos is not initialized" << endl;
    success = false;
  }

#ifdef HAVE_TPETRACORE_MPI
  // This label only gets used in an MPI build.
  // Some compilers warn on unused labels.
 EndOfTest:
#endif // HAVE_TPETRACORE_MPI

  // Clean up, if possible.
#ifdef HAVE_TPETRACORE_MPI
  if (mpiInitialized () && ! mpiFinalized ()) {
    (void) MPI_Finalize ();
  }
#endif // HAVE_TPETRACORE_MPI
  if (kokkosInitialized ()) {
    Kokkos::finalize ();
  }

  if (success) {
    std::cout << "End Result: TEST PASSED" << endl;
    return EXIT_SUCCESS;
  }
  else {
    std::cout << "End Result: TEST FAILED" << endl;
    return EXIT_FAILURE;
  }
}

