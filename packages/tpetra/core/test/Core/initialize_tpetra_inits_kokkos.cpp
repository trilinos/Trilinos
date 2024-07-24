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

void testMain (bool& success, int argc, char* argv[])
{
  using std::cout;
  using std::endl;

  // In this example, Tpetra::initialize is responsible for calling
  // Kokkos::initialize and Kokkos::finalize.
  if (Kokkos::is_initialized ()) {
    success = false;
    cout << "Kokkos::is_initialized() is true, "
      "before Tpetra::initialize was called." << endl;
    return;
  }
  Tpetra::initialize (&argc, &argv);

  if (! Kokkos::is_initialized ()) {
    success = false;
    cout << "Kokkos::is_initialized() is false, "
      "after Tpetra::initialize was called." << endl;
  }
  if (! Tpetra::isInitialized ()) {
    success = false;
    cout << "Tpetra::isInitialized() is false, "
      "even after Tpetra::initialize was called." << endl;
  }

  auto comm = Tpetra::getDefaultComm ();
  if (comm.is_null ()) {
    success = false;
    cout << "Tpetra::getDefaultComm() is null." << endl;
  }

  cout << "About to call Tpetra::finalize." << endl;
  Tpetra::finalize ();
  cout << "Called Tpetra::finalize." << endl;

  // Kokkos is like Tpetra; Kokkos::is_initialized() means "was
  // initialized and was not finalized."  That differs from MPI, where
  // MPI_Initialized only refers to MPI_Init and MPI_Finalized only
  // refers to MPI_Finalize.
  if (Kokkos::is_initialized ()) {
    success = false;
    cout << "Tpetra::finalize did not call Kokkos::finalize." << endl;
    return;
  }

  // MPI is no longer initialized, so we can't all-reduce on this.
  if (Tpetra::isInitialized ()) {
    success = false;
    cout << "Tpetra::isInitialized() is true, "
      "even after Tpetra::finalize has been called" << endl;
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
