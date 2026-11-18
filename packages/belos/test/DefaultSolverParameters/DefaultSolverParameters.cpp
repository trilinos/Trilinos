// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "BelosTypes.hpp"

#include <cmath>

// On Windows shared builds these static data members live in belos.dll and are
// dllimport'd here, so this test fails to link if belos does not export them.

namespace {

TEUCHOS_UNIT_TEST( DefaultSolverParameters, LinkAndValues )
{
  using Belos::DefaultSolverParameters;

  const double* const params[] = {
    &DefaultSolverParameters::convTol,
    &DefaultSolverParameters::polyTol,
    &DefaultSolverParameters::orthoKappa,
    &DefaultSolverParameters::resScaleFactor,
    &DefaultSolverParameters::impTolScale
  };
  for (const double* p : params) {
    TEST_ASSERT( p != nullptr );
    TEST_ASSERT( std::isfinite (*p) );
  }

  TEST_ASSERT( DefaultSolverParameters::convTol > 0.0 );
  TEST_ASSERT( DefaultSolverParameters::polyTol > 0.0 );
  TEST_ASSERT( DefaultSolverParameters::resScaleFactor > 0.0 );
  TEST_ASSERT( DefaultSolverParameters::impTolScale > 0.0 );
}

} // namespace
