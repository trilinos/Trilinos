// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Map.hpp"
#include "Kokkos_Core.hpp"
#include <iostream>
#include <cstdlib> // EXIT_SUCCESS, EXIT_FAILURE

// Test that Map does not try to initialize Kokkos again if Kokkos has
// already been initialized, and that Map does not attempt to
// double-finalize Kokkos in this case.

int
main (int argc, char* argv[])
{
  using std::cerr;
  using std::endl;

  if (Kokkos::is_initialized ()) {
    cerr << "FAILED: Before calling Kokkos::initialize, "
      "Kokkos reports that it has been initialized!" << endl;
    return EXIT_FAILURE;
  }
  Kokkos::initialize (argc, argv);
  if (! Kokkos::is_initialized ()) {
    cerr << "FAILED: After calling Kokkos::initialize, "
      "Kokkos is still not initialized!" << endl;
    return EXIT_FAILURE;
  }

  // Tpetra objects must never be created at main() scope, since they
  // may contain MPI and/or Kokkos objects that all must be freed
  // before MPI_Finalize resp. Kokkos::finalize have been called.
  {
    // Create a Map instance.  Make sure that this does not make Kokkos
    // raise an exception, e.g., due to double initialization.
    Tpetra::Map<> map1;

    if (! Kokkos::is_initialized ()) {
      cerr << "FAILED: After calling Kokkos::initialize, "
        "and after calling Map's constructor once, "
        "Kokkos is still not initialized!" << endl;
      return EXIT_FAILURE;
    }
    {
      Tpetra::Map<> map2; // inner scope, so its destructor gets invoked
      if (! Kokkos::is_initialized ()) {
        cerr << "FAILED: After calling Kokkos::initialize, "
          "and after calling Map's constructor twice, "
          "Kokkos is still not initialized!" << endl;
        return EXIT_FAILURE;
      }
    }
    if (! Kokkos::is_initialized ()) {
      cerr << "FAILED: After calling Kokkos::initialize, "
        "after calling Map's constructor twice, "
        "and after calling Map's destructor once, "
        "Kokkos is still not initialized!" << endl;
      return EXIT_FAILURE;
    }
  }

  Kokkos::finalize ();
  // Don't print PASSED!  That will make the test think that it
  // passed, before any atexit hooks may have been called (they should
  // not be called in this case, but we want to test that they haven't
  // been called).  We use FAIL_REGULAR_EXPRESSION (see
  // CMakeLists.txt) to help with this.
  return EXIT_SUCCESS;
}

