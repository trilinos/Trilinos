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
#include <cstdlib> // std::getenv

namespace {

/*
 * Note: The tests for the Behavior class are scattered across several files,
 * rather than being confined to a single file.  The reason is that the Behavior
 * class is instantiated one time only and once environment variables are read,
 * they are cached for future use.  Therefore, to test several values of an
 * environment variable, several tests need to be created (one for each distinct
 * value of the environment variable).
*/

TEUCHOS_UNIT_TEST(Behavior, Default)
{
  bool verbose_default = false;
  bool verb = Tpetra::Details::Behavior::verbose();
  TEUCHOS_TEST_ASSERT(verb==verbose_default, out, success);

  // Print current values of other behaviors. 
  // Can't test against default since these behaviors may be 
  // changed by environment variables (in which case, test against
  // default fails)
  std::cout << "\n        GPU-aware MPI?  " 
            << Tpetra::Details::Behavior::assumeMpiIsGPUAware()
            << "\n";

  std::cout << "\n        Tpetra Debug?  "
            << Tpetra::Details::Behavior::debug()
            << "\n";
}

TEUCHOS_UNIT_TEST(Behavior, verbosePrintCountThreshold) {
  // We only require that the default be between these values.
  const size_t maxVal (1000);
  const size_t minVal (100);
  const size_t val0 =
    Tpetra::Details::Behavior::verbosePrintCountThreshold();
  TEST_ASSERT( val0 >= minVal && val0 <= maxVal );

  const size_t val1 =
    Tpetra::Details::Behavior::verbosePrintCountThreshold();
  TEST_ASSERT( val1 >= minVal && val1 <= maxVal );
}

} // namespace (anonymous)
