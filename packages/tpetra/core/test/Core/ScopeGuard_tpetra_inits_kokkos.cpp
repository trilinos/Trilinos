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

  // In this example, Tpetra::ScopeGuard is responsible for calling
  // Kokkos::initialize and Kokkos::finalize.
  if (Kokkos::is_initialized ()) {
    success = false;
    cout << "Kokkos::is_initialized() is true, "
      "before Tpetra::ScopeGuard was created." << endl;
    return;
  }
  {
    Tpetra::ScopeGuard tpetraScope (&argc, &argv);

    if (! Kokkos::is_initialized ()) {
      success = false;
      cout << "Kokkos::is_initialized() is false, "
        "after Tpetra::ScopeGuard was created." << endl;
    }
    if (! Tpetra::isInitialized ()) {
      success = false;
      cout << "Tpetra::isInitialized() is false, "
        "even after Tpetra::ScopeGuard was created." << endl;
    }

    auto comm = Tpetra::getDefaultComm ();
    if (comm.is_null ()) {
      success = false;
      cout << "Tpetra::getDefaultComm() is null." << endl;
    }

    cout << "About to leave Tpetra scope." << endl;
  }
  cout << "Left Tpetra scope." << endl;

  // Kokkos is like Tpetra; Kokkos::is_initialized() means "was
  // initialized and was not finalized."  That differs from MPI, where
  // MPI_Initialized only refers to MPI_Init and MPI_Finalized only
  // refers to MPI_Finalize.
  if (Kokkos::is_initialized ()) {
    success = false;
    cout << "Tpetra::ScopeGuard::~ScopeGuard did not call Kokkos::finalize." << endl;
    return;
  }

  // MPI is no longer initialized, so we can't all-reduce on this.
  if (Tpetra::isInitialized ()) {
    success = false;
    cout << "Tpetra::isInitialized() is true, "
      "even after Tpetra::ScopeGuard::~ScopeGuard has been called" << endl;
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
