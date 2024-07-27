// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>
#include <string>
#include <cstdlib>
#include <Tpetra_Core.hpp>
#include <Tpetra_TestingUtilities.hpp>
#include <Tpetra_Details_Behavior.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_CommHelpers.hpp>

/*
 * Note: The tests for the Behavior class are scattered across several files,
 * rather than being confined to a single file.  The reason is that the Behavior
 * class is instantiated one time only and once environment variables are read,
 * they are cached for future use.  Therefore, to test several values of an
 * environment variable, several tests need to be created (one for each distinct
 * value of the environment variable).
*/

namespace {

TEUCHOS_STATIC_SETUP()
{
  setenv("TPETRA_DEBUG", "ON", 1);
  setenv("TPETRA_VERBOSE", "ON", 1);
  setenv("TPETRA_ASSUME_GPU_AWARE_MPI", "ON", 1);
  setenv("CUDA_LAUNCH_BLOCKING", "1", 1);
}

TEUCHOS_UNIT_TEST(Behavior, On)
{

  // TPETRA_DEBUG was set globally in TEUCHOS_STATIC_SETUP to ON, so any query
  // on TPETRA_DEBUG should evaluate to true, including named variants.
  bool dbg = Tpetra::Details::Behavior::debug();
  TEUCHOS_TEST_ASSERT(dbg, out, success);
  bool dbg_named = Tpetra::Details::Behavior::debug("Named");
  TEUCHOS_TEST_ASSERT(dbg_named, out, success);

  // TPETRA_VERBOSE was set globally in TEUCHOS_STATIC_SETUP to ON, so any query
  // on TPETRA_VERBOSE should evaluate to true, including named variants.
  bool verb = Tpetra::Details::Behavior::verbose();
  TEUCHOS_TEST_ASSERT(verb, out, success);
  bool verb_named = Tpetra::Details::Behavior::verbose("Named");
  TEUCHOS_TEST_ASSERT(verb_named, out, success);

  // TPETRA_ASSUME_GPU_AWARE_MPI was set globally in TEUCHOS_STATIC_SETUP to ON,
  // so any query on TPETRA_ASSUME_GPU_AWARE_MPI should evaluate to true
  bool gpu_aware_mpi = Tpetra::Details::Behavior::assumeMpiIsGPUAware();
  TEUCHOS_TEST_ASSERT(gpu_aware_mpi, out, success);

  TEST_ASSERT(Tpetra::Details::Behavior::cudaLaunchBlocking());
}
} // namespace (anonymous)
