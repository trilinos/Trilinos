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

// Test that Map's constructor initializes Kokkos if it is not already
// initialized.  It must also finalize Kokkos at exit, but the only
// portable way to test that would be to run Valgrind and ensure no
// memory leaks.

int
main (int argc, char* argv[])
{
  using std::cerr;
  using std::endl;

  if (Kokkos::is_initialized ()) {
    cerr << "FAILED: Before calling Map's constructor, "
      "and without having called Kokkos::initialize ourselves, "
      "Kokkos reports that it has been initialized!" << endl;
    return EXIT_FAILURE;
  }

  {
    // Create a Map instance.  Make sure that this does not make
    // Kokkos raise an exception, e.g., due to double initialization.
    // Do this in an inner (within main()) scope, so we can test that
    // case.  Another test in this directory exercises the case where
    // node instances are always at main() scope.
    Tpetra::Map<> map1;

    if (! Kokkos::is_initialized ()) {
      cerr << "FAILED: After calling Kokkos::initialize, "
        "and after calling Map's constructor once, "
        "Kokkos is still not initialized!" << endl;
      return EXIT_FAILURE;
    }

    // Kokkos itself will check that the execution space has been
    // initialized, if we create a View, since Kokkos will run a
    // parallel_for to initialize the View.
    try {
      using device_type = Tpetra::Map<>::device_type;
      Kokkos::View<double*, device_type> testView ("testView", 10);
    }
    catch (...) {
      cerr << "FAILED: After calling Kokkos::initialize, "
        "and after calling Map's constructor once, "
        "Kokkos::View creation fails, "
        "likely because Kokkos is still not initialized!" << endl;
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

  // Don't print PASSED!  That will make the test think that it
  // passed, even though things may still be happening (due to atexit)
  // after main() exits.  We use FAIL_REGULAR_EXPRESSION (see
  // CMakeLists.txt) so that the test fails if and only if it prints
  // "FAILED:".
  return EXIT_SUCCESS;
}

