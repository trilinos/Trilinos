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

  if (Kokkos::is_initialized ()) {
    success = false;
    cout << "Kokkos::is_initialized() is true, "
      "even before Kokkos::initialize was called." << endl;
    return;
  }

  Teuchos::GlobalMPISession mpiSession(&argc, &argv); // before Kokkos::initialize
  Kokkos::initialize (argc, argv);
  if (! Kokkos::is_initialized ()) {
    success = false;
    cout << "Kokkos::is_initialized() is false, "
      "even after Kokkos::initialize was called." << endl;
    return;
  }

  // In this example, the "user" has called Kokkos::initialize before
  // Tpetra::initialize is called.  Tpetra::initialize must not try to
  // call it again.
  Tpetra::initialize (&argc, &argv);
  if (! Tpetra::isInitialized ()) {
    success = false;
    cout << "Tpetra::isInitialized() is false, "
      "even after Tpetra::initialize was called."
      << endl;
    return;
  }
  if (! Kokkos::is_initialized ()) {
    success = false;
    cout << "Kokkos::is_initialized() is false, "
      "even after Kokkos::initialize and Tpetra::initialize were called."
      << endl;
    // Let the program keep going, so MPI (if applicable) gets finalized.
  }

  cout << "About to call Tpetra::finalize" << endl;
  Tpetra::finalize ();
  cout << "Called Tpetra::finalize" << endl;
  if (Tpetra::isInitialized ()) {
    success = false;
    cout << "Tpetra::isInitialized() is true, "
      "even after Tpetra::finalize was called."
      << endl;
    // Let the program keep going, so Kokkos gets finalized.
  }

  // Since the "user" is responsible for calling Kokkos::finalize,
  // Tpetra::finalize should NOT have called Kokkos::finalize.
  if (! Kokkos::is_initialized ()) {
    success = false;
    cout << "Kokkos::is_initialized() is false, "
      "after Tpetra::initialize was called." << endl;
  }
  Kokkos::finalize ();
  if (Kokkos::is_initialized ()) {
    success = false;
    cout << "Kokkos::is_initialized() is true, "
      "even after Kokkos::initialize was called." << endl;
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
